#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HclusterParseOutput

=cut

=head1 SYNOPSIS

my $aa = $sdba->get_AnalysisAdaptor;
my $analysis = $aa->fetch_by_logic_name('HclusterParseOutput');
my $rdb = new Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HclusterParseOutput(
                         -input_id   => "{'mlss_id'=>40069}",
                         -analysis   => $analysis);

$rdb->fetch_input
$rdb->run;

=cut

=head1 DESCRIPTION

This is the RunnableDB that parses the output of Hcluster, stores the clusters as trees without internal structure
(each tree will have one root and several leaves) and dataflows the cluster_ids down branch #2.

=cut

=head1 CONTACT

  Contact Albert Vilella on module implemetation/design detail: avilella@ebi.ac.uk
  Contact Abel Ureta-Vidal on EnsEMBL/Compara: abel@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::HclusterParseOutput;

use strict;
use Bio::EnsEMBL::Compara::NestedSet;
use Bio::EnsEMBL::Compara::Graph::ConnectedComponents;
use Bio::EnsEMBL::Compara::MethodLinkSpeciesSet;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');


sub param_defaults {
    return {
        'clusterset_id'         => 1,
    };
}


sub fetch_input {
  my $self = shift;

  my $mlss_id   = $self->param('mlss_id') or die "'mlss_id' is an obligatory parameter";

  my $protein_tree_adaptor = $self->compara_dba->get_ProteinTreeAdaptor;
  my $member_adaptor       = $self->compara_dba->get_MemberAdaptor;

  my $cluster_dir   = $self->param('cluster_dir');
  my $filename      = $cluster_dir . '/hcluster.out';

  # FIXME: load the entire file in a hash and store in decreasing
  # order by cluster size this will make big clusters go first in the
  # alignment process, which makes sense since they are going to take
  # longer to process anyway
  my $clusterset;
  $clusterset = $protein_tree_adaptor->fetch_node_by_node_id($self->param('clusterset_id'));
  if (!defined($clusterset)) {
    $self->param('ccEngine', Bio::EnsEMBL::Compara::Graph::ConnectedComponents->new() );
    $clusterset = $self->param('ccEngine')->clusterset;
    $self->throw("no clusters generated") unless($clusterset);

    $clusterset->name("PROTEIN_TREES");
    $protein_tree_adaptor->store_node($clusterset);
    printf("clusterset_id %d\n", $clusterset->node_id);
    $self->param('clusterset_id', $clusterset->node_id);
  }

  open(FILE, $filename) or die "Could not open '$filename' for reading : $!";
  while (<FILE>) {
    # 0       0       0       1.000   2       1       697136_68,
    # 1       0       39      1.000   3       5       1213317_31,1135561_22,288182_42,426893_62,941130_38,
    chomp $_;

    my ($cluster_id, $dummy1, $dummy2, $dummy3, $dummy4, $dummy5, $cluster_list) = split("\t",$_);

    next if ($dummy5 < 2);
    $cluster_list =~ s/\,^//;
    my @cluster_list = split(",",$cluster_list);

    # If it's a singleton, we don't store it as a protein tree
    next if (2 > scalar(@cluster_list));

    my $cluster = new Bio::EnsEMBL::Compara::NestedSet;
    $clusterset->add_child($cluster);

    foreach my $member_hcluster_id (@cluster_list) {
      my ($pmember_id,$genome_db_id) = split("_",$member_hcluster_id);

      my $node = new Bio::EnsEMBL::Compara::NestedSet;
      $node->node_id($pmember_id);
      $cluster->add_child($node);
      $cluster->clusterset_id($self->param('clusterset_id'));
      #leaves are NestedSet objects, bless to make into AlignedMember objects
      bless $node, "Bio::EnsEMBL::Compara::AlignedMember";

      #the building method uses member_id's to reference unique nodes
      #which are stored in the node_id value, copy to member_id
      $node->member_id($node->node_id);
      $node->method_link_species_set_id($mlss_id);
    }

        # Store the cluster:
    $protein_tree_adaptor->store($cluster);

    my $leafcount = scalar(@{$cluster->get_all_leaves});
    $cluster->store_tag('gene_count', $leafcount);

  }
  close FILE;
}


sub write_output {
    my $self = shift;

    my $clusterset = $self->compara_dba->get_ProteinTreeAdaptor->fetch_node_by_node_id($self->param('clusterset_id'));
    if (!defined($clusterset)) {
        $clusterset = $self->param('ccEngine')->clusterset;
    }

    foreach my $cluster (@{$clusterset->children()}) {
        $self->dataflow_output_id({
            'protein_tree_id'   => $cluster->node_id(),
            'clusterset_id'     => $clusterset->node_id(),
        }, 2);
    }
}

1;
