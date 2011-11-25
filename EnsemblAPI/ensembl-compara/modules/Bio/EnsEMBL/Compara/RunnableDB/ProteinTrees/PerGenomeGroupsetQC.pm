#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::PerGenomeGroupsetQC

=cut

=head1 SYNOPSIS

my $db           = Bio::EnsEMBL::Compara::DBAdaptor->new($locator);
my $sillytemplate = Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::PerGenomeGroupsetQC->new
  (
   -db         => $db,
   -input_id   => $input_id,
   -analysis   => $analysis
  );
$sillytemplate->fetch_input(); #reads from DB
$sillytemplate->run();
$sillytemplate->output();
$sillytemplate->write_output(); #writes to DB

=cut


=head1 DESCRIPTION

This Analysis will take the sequences from a cluster, the cm from
nc_profile and run a profiled alignment, storing the results as
cigar_lines for each sequence.

=cut


=head1 CONTACT

  Contact Albert Vilella on module implementation/design detail: avilella@ebi.ac.uk
  Contact Ewan Birney on EnsEMBL in general: birney@sanger.ac.uk

=cut


=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut


package Bio::EnsEMBL::Compara::RunnableDB::ProteinTrees::PerGenomeGroupsetQC;

use strict;

use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');


sub run {
    my $self = shift @_;

    my $genome_db_id            = $self->param('genome_db_id') or die "'genome_db_id' is an obligatory parameter";
    my $groupset_tag            = $self->param('groupset_tag') or die "'groupset_tag' is an obligatory parameter";
    my $groupset_node           = $self->compara_dba->get_ProteinTreeAdaptor->fetch_all_roots->[0] or die "Could not fetch groupset node";

    my $this_orphans            = $self->fetch_gdb_orphan_genes($self->compara_dba, $genome_db_id);
    my $total_orphans_num       = scalar keys (%$this_orphans);
    my $total_num_genes         = scalar @{ $self->compara_dba->get_MemberAdaptor->fetch_all_by_source_genome_db_id('ENSEMBLGENE',$genome_db_id) };
    my $proportion_orphan_genes = $total_orphans_num/$total_num_genes;

    $groupset_node->store_tag("$genome_db_id".'_total_orphans_num' . '_' . $groupset_tag, $total_orphans_num);
    $groupset_node->store_tag("$genome_db_id".'_prop_orphans' . '_' . $groupset_tag, $proportion_orphan_genes);

    return unless $self->param('reuse_this');

    my $reuse_db                = $self->param('reuse_db') or die "'reuse_db' connection parameters hash has to be defined in reuse mode";

    my $reuse_compara_dba       = $self->go_figure_compara_dba($reuse_db);    # may die if bad parameters

    my $reuse_orphans           = $self->fetch_gdb_orphan_genes($reuse_compara_dba, $genome_db_id);
    my %common_orphans = ();
    my %new_orphans = ();
    foreach my $this_orphan_id (keys %$this_orphans) {
        if($reuse_orphans->{$this_orphan_id}) {
            $common_orphans{$this_orphan_id} = 1;
        } else {
            $new_orphans{$this_orphan_id} = 1;
        }
    }
    my $common_orphans_num = scalar keys (%common_orphans);
    my $new_orphans_num    = scalar keys (%new_orphans);

    $groupset_node->store_tag("$genome_db_id".'_common_orphans_num' . '_' . $groupset_tag, $common_orphans_num);
    $groupset_node->store_tag("$genome_db_id".'_new_orphans_num' . '_' . $groupset_tag, $new_orphans_num);
}


sub fetch_gdb_orphan_genes {
    my ($self, $given_compara_dba, $genome_db_id) = @_;

    my %orphan_stable_id_hash = ();

    my $sql = "SELECT m3.stable_id from member m2, member m3, subset_member sm where m3.member_id=m2.gene_member_id and m2.source_name='ENSEMBLPEP' and sm.member_id=m2.member_id and sm.member_id in (SELECT m1.member_id from member m1 left join protein_tree_member ptm on m1.member_id=ptm.member_id where ptm.member_id IS NULL and m1.genome_db_id=$genome_db_id)";

    my $sth = $given_compara_dba->dbc->prepare($sql);
    $sth->execute();

    while(my ($member) = $sth->fetchrow()) {
        $orphan_stable_id_hash{$member} = 1;
    }

    return \%orphan_stable_id_hash;
}

1;
