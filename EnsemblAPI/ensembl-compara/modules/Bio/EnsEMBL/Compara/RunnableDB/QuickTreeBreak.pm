#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::QuickTreeBreak

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $quicktreebreak = Bio::EnsEMBL::Compara::RunnableDB::QuickTreeBreak->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$quicktreebreak->fetch_input(); #reads from DB
$quicktreebreak->run();
$quicktreebreak->output();
$quicktreebreak->write_output(); #writes to DB

=cut


=head1 DESCRIPTION

This Analysis/RunnableDB is designed to take ProteinTree as input.

This must already have a multiple alignment run on it. It uses that
alignment as input into the QuickTree program which then generates a
simple phylogenetic tree to be broken down into 2 pieces.

Google QuickTree to get the latest tar.gz from the Sanger.
Google sreformat to get the sequence reformatter that switches from fasta to stockholm.

input_id/parameters format eg: "{'protein_tree_id'=>1234,'clusterset_id'=>1}"
    protein_tree_id : use 'id' to fetch a cluster from the ProteinTree

=cut


=head1 CONTACT

  Contact Albert Vilella on module implementation/design detail: avilella@ebi.ac.uk
  Contact Javier Herrero on EnsEMBL/Compara: jherrero@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Compara::RunnableDB::QuickTreeBreak;

use strict;
use IO::File;
use File::Basename;
use Time::HiRes qw(time gettimeofday tv_interval);

use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::Graph::NewickParser;
use Bio::SimpleAlign;
use Bio::AlignIO;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');


sub param_defaults {
    return {
        'sreformat_exe'         => '/usr/local/ensembl/bin/sreformat',
    };
}


=head2 fetch_input

    Title   :   fetch_input
    Usage   :   $self->fetch_input
    Function:   Fetches input data for repeatmasker from the database
    Returns :   none
    Args    :   none

=cut


sub fetch_input {
    my $self = shift @_;

    $self->check_if_exit_cleanly;

    my $protein_tree_id     = $self->param('protein_tree_id') or die "'protein_tree_id' is an obligatory parameter";
    my $protein_tree        = $self->compara_dba->get_ProteinTreeAdaptor->fetch_node_by_node_id( $protein_tree_id )
                                        or die "Could not fetch protein_tree with protein_tree_id='$protein_tree_id'";
    $self->param('protein_tree', $protein_tree);


    # FIXME: This can be quite dangerous, as there may be more than one in the DB by accident:
      #
      # We fetch the mlssID that is later needed for the newly stored leaves
    my @protein_trees_mlsses = @{$self->compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_all_by_method_link_type('PROTEIN_TREES')};
    $self->param('mlss_id', $protein_trees_mlsses[0]->dbID || 0 );
}


=head2 run

    Title   :   run
    Usage   :   $self->run
    Function:   runs NJTREE PHYML
    Returns :   none
    Args    :   none

=cut


sub run {
    my $self = shift @_;

    $self->check_if_exit_cleanly;
    $self->run_quicktreebreak;
}


=head2 write_output

    Title   :   write_output
    Usage   :   $self->write_output
    Function:   stores proteintree
    Returns :   none
    Args    :   none

=cut


sub write_output {
    my $self = shift @_;

    $self->check_if_exit_cleanly;
    $self->store_proteintrees;
}


sub DESTROY {
  my $self = shift;

  if($self->param('protein_tree')) {
    printf("QuickTreeBreak::DESTROY releasing tree\n") if($self->debug);

    $self->param('protein_tree')->release_tree;
    $self->param('protein_tree', undef);

    $self->param('max_subtree')->release_tree;
    $self->param('new_subtree')->release_tree;
    $self->param('remaining_subtree')->release_tree;

    $self->param('max_subtree', undef);
    $self->param('new_subtree', undef);
    $self->param('remaining_subtree', undef);
  }

  $self->SUPER::DESTROY if $self->can("SUPER::DESTROY");
}


##########################################
#
# internal methods
#
##########################################


sub run_quicktreebreak {
  my $self = shift;

  my $starttime = time()*1000;

  my $input_aln = $self->dumpTreeMultipleAlignmentToWorkdir ( $self->param('protein_tree') ) or return;

  my $quicktreebreak_executable = $self->analysis->program_file || '';

  unless (-e $quicktreebreak_executable) {
    $quicktreebreak_executable = "/nfs/users/nfs_a/avilella/src/quicktree/quicktree_1.1/bin/quicktree";
  }

  $self->throw("can't find a quicktree executable to run. Tried $quicktreebreak_executable \n") 
    unless(-e $quicktreebreak_executable);

  my $cmd = $quicktreebreak_executable;
  $cmd .= " -out t -in a";
  $cmd .= " ". $input_aln;

  #/nfs/users/nfs_a/avilella/src/quicktree/quicktree_1.1/bin/quicktree -out t
  # -in a /tmp/worker.12270/proteintree_517373.stk

  $self->compara_dba->dbc->disconnect_when_inactive(1);
  print("$cmd\n") if($self->debug);
  open(RUN, "$cmd |") or $self->throw("error running quicktree, $!\n");
  my @output = <RUN>;
  my $exit_status = close(RUN);
  if (!$exit_status) {
    $self->throw("error running quicktree, $!\n");
  }
  $self->compara_dba->dbc->disconnect_when_inactive(0);

  my $quicktree_newick_string = '';
  foreach my $line (@output) {
    $line =~ s/\n//;
    $quicktree_newick_string .= $line;
  }

  #parse the tree into the datastucture
  $self->generate_subtrees( $quicktree_newick_string );

  my $runtime = time()*1000-$starttime;
  $self->param('protein_tree')->store_tag('QuickTreeBreak_runtime_msec', $runtime);
}


########################################################
#
# ProteinTree input/output section
#
########################################################

sub dumpTreeMultipleAlignmentToWorkdir {
  my $self = shift;
  my $protein_tree = shift;

  $self->param('original_leafcount', scalar(@{$protein_tree->get_all_leaves}) );
  if($self->param('original_leafcount')<3) {
    printf(STDERR "tree cluster %d has <3 proteins - can not build a tree\n", $protein_tree->node_id);
    return undef;
  }

  my $file_root = $self->worker_temp_directory. "proteintree_". $protein_tree->node_id;
  $file_root =~ s/\/\//\//g;  # converts any // in path to /

  my $aln_file = $file_root . '.aln';
  return $aln_file if(-e $aln_file);
  if($self->debug) {
    printf("dumpTreeMultipleAlignmentToWorkdir : %d members\n", $self->param('original_leafcount') );
    print("aln_file = '$aln_file'\n");
  }

  open(OUTSEQ, ">$aln_file") or $self->throw("Error opening $aln_file for write");

  # Using append_taxon_id will give nice seqnames_taxonids needed for
  # njtree species_tree matching
  my %sa_params = ($self->param('use_genomedb_id')) ? ('-APPEND_GENOMEDB_ID', 1) : ('-APPEND_TAXON_ID', 1);

  my $sa = $protein_tree->get_SimpleAlign ( -id_type => 'MEMBER', %sa_params );
  $sa->set_displayname_flat(1);
  my $alignIO = Bio::AlignIO->newFh ( -fh => \*OUTSEQ, -format => "fasta" );
  print $alignIO $sa;

  close OUTSEQ;

  print STDERR "Using sreformat to change to stockholm format\n" if ($self->debug);
  my $stk_file = $file_root . '.stk';
  
  my $sreformat_exe = $self->param('sreformat_exe');
  
  my $cmd = "$sreformat_exe stockholm $aln_file > $stk_file";

  unless( system("$cmd") == 0) {
    print("$cmd\n");
    $self->throw("error running sreformat with cmd $cmd: $!\n");
  }

  return $stk_file;
}


sub store_proteintrees {
  my $self = shift;

  $self->delete_original_cluster;
  $self->store_clusters;

  if($self->debug >1) {
    print("done storing\n");
  }

  return undef;
}

sub store_clusters {
  my $self = shift;

  my $protein_tree_adaptor = $self->compara_dba->get_ProteinTreeAdaptor;
  my $starttime = time();

  my $clusterset = $protein_tree_adaptor->fetch_node_by_node_id($self->param('clusterset_id'));
  $self->throw("no clusterset found: $!\n") unless($clusterset);

  $clusterset->no_autoload_children; # Important so that only the two below are used
  $clusterset->add_child($self->param('new_subtree'));
  $clusterset->add_child($self->param('remaining_subtree'));

  my $clusters = $clusterset->children;
  foreach my $cluster (@{$clusters}) {
    my $node_id = $protein_tree_adaptor->store($cluster);
    # Although the leaves wont have the right root_id pointing to the $cluster->node_id,
    # this will be solved when we store back the results after the new MSA job.

    #calc residue count total
    my $leafcount = scalar(@{$cluster->get_all_leaves});
    $cluster->store_tag('gene_count', $leafcount);
    $cluster->store_tag('original_cluster', $self->param('original_cluster')->node_id);
    print STDERR "Stored $node_id with $leafcount leaves\n" if ($self->debug);

    next if($leafcount<2);

    # Dataflow clusters
    # This will create a new MSA alignment job for each of the newly generated clusters
    my $output_id = sprintf("{'protein_tree_id'=>%d, 'clusterset_id'=>%d}", 
                            $node_id, $clusterset->node_id);

    $self->dataflow_output_id($output_id, 1);
    print STDERR "Created new cluster $node_id\n";
  }
}

sub delete_original_cluster {
  my $self = shift;

  my $delete_time = time();
  my $original_cluster_node_id = $self->param('original_cluster')->node_id;
  $self->delete_old_orthotree_tags;

  my $tree_node_id = $original_cluster_node_id;
  my $sql1 = "delete h.*, hm.* from homology h, homology_member hm where h.homology_id=hm.homology_id and h.tree_node_id=$tree_node_id";
  my $sth1 = $self->dbc->prepare($sql1);
  $sth1->execute;
  $sth1->finish;

  $self->param('original_cluster')->adaptor->store_supertree_node_and_under( $self->param('original_cluster') );
  printf("%1.3f secs to copy old cluster $original_cluster_node_id into supertree tables\n", time()-$delete_time);
  $self->param('original_cluster')->adaptor->delete_node_and_under($self->param('original_cluster'));
  printf("%1.3f secs to delete old cluster $original_cluster_node_id\n", time()-$delete_time);

  return 1;

}

sub delete_old_orthotree_tags {
  my $self = shift;

  print "deleting old orthotree tags\n" if ($self->debug);
  my @node_ids;
  my $left_index  = $self->param('protein_tree')->left_index;
  my $right_index = $self->param('protein_tree')->right_index;
  my $tree_root_node_id = $self->param('protein_tree')->node_id;
  # Include the root_id as well as the rest of the nodes within the tree
  push @node_ids, $tree_root_node_id;
  my $sql = "select ptn.node_id from protein_tree_node ptn where ptn.left_index>$left_index and ptn.right_index<$right_index";
  my $sth = $self->dbc->prepare($sql);
  $sth->execute;
  while (my $aref = $sth->fetchrow_arrayref) {
    my ($node_id) = @$aref;
    push @node_ids, $node_id;
  }

  my @list_ids;
  foreach my $id (@node_ids) {
    push @list_ids, $id;
    if (scalar @list_ids == 2000) {
      my $sql = "delete from protein_tree_tag where node_id in (".join(",",@list_ids).") and tag in ('duplication_confidence_score','taxon_id','taxon_name','OrthoTree_runtime_msec','OrthoTree_types_hashstr')";
      my $sth = $self->dbc->prepare($sql);
      $sth->execute;
      $sth->finish;
      @list_ids = ();
    }
  }
  
  if (scalar @list_ids) {
    my $sql = "delete from protein_tree_tag where node_id in (".join(",",@list_ids).") and tag in ('duplication_confidence_score','taxon_id','taxon_name','OrthoTree_runtime_msec','OrthoTree_types_hashstr')";
    my $sth = $self->dbc->prepare($sql);
    $sth->execute;
    $sth->finish;
    @list_ids = ();
  }

  return undef;
}

sub generate_subtrees {
    my $self                    = shift @_;
    my $quicktree_newick_string = shift @_;

    my $mlss_id = $self->param('mlss_id');
    my $protein_tree = $self->param('protein_tree');

  #cleanup old tree structure- 
  #  flatten and reduce to only AlignedMember leaves
  $protein_tree->flatten_tree;
  $protein_tree->print_tree(20) if($self->debug);
  foreach my $node (@{$protein_tree->get_all_leaves}) {
    next if($node->isa('Bio::EnsEMBL::Compara::AlignedMember'));
    $node->disavow_parent;
  }

  #parse newick into a new tree object structure
  my $newtree = Bio::EnsEMBL::Compara::Graph::NewickParser::parse_newick_into_tree($quicktree_newick_string);
  $newtree->print_tree(20) if($self->debug > 1);
  # get rid of the taxon_id needed by njtree -- name tag
  foreach my $leaf (@{$newtree->get_all_leaves}) {
    my $quicktreebreak_name = $leaf->get_tagvalue('name');
    $quicktreebreak_name =~ /(\d+)\_\d+/;
    my $member_name = $1;
    $leaf->add_tag('name', $member_name);
    bless $leaf, "Bio::EnsEMBL::Compara::AlignedMember";
    $leaf->member_id($member_name);
  }

  # Break the tree by immediate children recursively
  my @children;
  my $keep_breaking = 1;
  $self->param('max_subtree', $newtree);
  while ($keep_breaking) {
    @children = @{$self->param('max_subtree')->children};
    my $max_num_leaves = 0;
    foreach my $child (@children) {
      my $num_leaves = scalar(@{$child->get_all_leaves});
      if ($num_leaves > $max_num_leaves) {
        $max_num_leaves = $num_leaves;
        $self->param('max_subtree', $child);
      }
    }
    # Broke down to half, happy with it
    my $proportion = ($max_num_leaves*100/$self->param('original_leafcount') );
    print STDERR "QuickTreeBreak iterate -- $max_num_leaves ($proportion)\n" if ($self->debug);
    if ($proportion <= 50) {
      $keep_breaking = 0;
    }
  }

  # Create a copy of what is not max_subtree
  $self->param('remaining_subtree', $self->param('protein_tree')->copy);
  $self->param('new_subtree',       $self->param('protein_tree')->copy);
  $self->param('new_subtree')->flatten_tree;
  $self->param('remaining_subtree')->flatten_tree;
  my $subtree_leaves;
  foreach my $leaf (@{$self->param('max_subtree')->get_all_leaves}) {
    $subtree_leaves->{$leaf->member_id} = 1;
  }
  foreach my $leaf (@{$self->param('new_subtree')->get_all_leaves}) {
    unless (defined $subtree_leaves->{$leaf->member_id}) {
      print $leaf->name," leaf disavowing parent\n" if $self->debug;
      $leaf->disavow_parent;
    }
    $leaf->method_link_species_set_id($self->param('mlss_id'));
  }
  foreach my $leaf (@{$self->param('remaining_subtree')->get_all_leaves}) {
    if (defined $subtree_leaves->{$leaf->member_id}) {
      print $leaf->name," leaf disavowing parent\n" if $self->debug;
      $leaf->disavow_parent;
    }
    $leaf->method_link_species_set_id($self->param('mlss_id'));
  }
  $self->param('remaining_subtree', $self->param('remaining_subtree')->minimize_tree);
  $self->param('new_subtree',       $self->param('new_subtree')->minimize_tree);

  # Some checks
  $self->throw("QuickTreeBreak: Failed to generate subtrees: $!\n")  unless(defined($self->param('new_subtree')) && defined($self->param('remaining_subtree')));
  my  $final_original_num = scalar @{$self->param('protein_tree')->get_all_leaves};
  my       $final_max_num = scalar @{$self->param('new_subtree')->get_all_leaves};
  my $final_remaining_num = scalar @{$self->param('remaining_subtree')->get_all_leaves};

  if(($final_max_num + $final_remaining_num) != $final_original_num) {
    $self->throw("QuickTreeBreak: Incorrect sum of leaves [$final_max_num + $final_remaining_num != $final_original_num]: $!\n");
  }

  $self->param('original_cluster', $self->param('protein_tree'));
  return undef;
}

1;
