package Bio::EnsEMBL::Compara::DBSQL::MemberAdaptor;

use strict; 
use warnings;
use Bio::EnsEMBL::Compara::Member;
use Bio::EnsEMBL::Compara::Attribute;
use Bio::EnsEMBL::Compara::DBSQL::SequenceAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning stack_trace_dump);

our @ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);



sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);

  $self->{'_member_cache'} = {};
  return $self;
}



=head2 list_internal_ids

  Arg        : None
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub list_internal_ids {
  my $self = shift;
  
  my ($name, $syn) = @{$self->tables->[0]};
  my $sql = "SELECT ${syn}.${name}_id from ${name} ${syn}";
  
  my $sth = $self->prepare($sql);
  $sth->execute;  
  
  my $internal_id;
  $sth->bind_columns(\$internal_id);

  my @internal_ids;
  while ($sth->fetch()) {
    push @internal_ids, $internal_id;
  }

  $sth->finish;

  return \@internal_ids;
}


sub member_cache {  
  my ($self,$id,$val) = @_; 

  my $result ; 
  if ( $id ) {  
    if ( defined $val ) {   
       $self->{_member_cache}{$id}=$val;
    }  
    $result =$self->{_member_cache}{$id};
  }   
  return $result ; 
}  

=head2 fetch_by_dbID

  Arg [1]    : int $id
               the unique database identifier for the feature to be obtained
  Example    : $feat = $adaptor->fetch_by_dbID(1234);
  Description: Returns the Member created from the database defined by the
               the id $id.
  Returntype : Bio::EnsEMBL::Compara::Member
  Exceptions : thrown if $id is not defined
  Caller     : general

=cut
sub fetch_by_dbID {
  my ($self,$id) = @_;

  unless(defined $id) {
    throw("fetch_by_dbID must have an id");
  }

## adapted function to cache the member 
  my $obj; 
  # check if member is cached already  
  my $member_is_cached = $self->member_cache($id) ; 

  if ( defined $member_is_cached ) {   
    $obj = $member_is_cached;
  } else { 
    my ($name, $syn) = @{$self->tables->[0]};
  
    #construct a constraint like 't1.table1_id = 1'
    my $constraint = "${syn}.${name}_id = $id";
  
    #return first element of _generic_fetch list
  
    ($obj) = @{$self->_generic_fetch($constraint)};  
    $self->member_cache($id,$obj);  
  }
  return $obj;
}



## previous function 
#   my ($name, $syn) = @{$self->tables->[0]};
#   #construct a constraint like 't1.table1_id = 1'
#   my $constraint = "${syn}.${name}_id = $id";
  #return first element of _generic_fetch list
#   my ($obj) = @{$self->_generic_fetch($constraint)};
#   return $obj;
#}






sub fetch_by_dbIDs {
  my $self = shift;

  my $ids = join(',' , @_);
  my $constraint = "m.member_id in ($ids)";
  return $self->_generic_fetch($constraint);
}

sub fetch_all_by_sequence_id {
  my ($self, $sequence_id) = @_;

  $self->throw("sequence_id arg is required\n")
    unless (defined($sequence_id));

  my $constraint = "m.sequence_id = $sequence_id";
  return $self->_generic_fetch($constraint);
}

=head2 fetch_by_source_stable_id

  Arg [1]    : (optional) string $source_name
  Arg [2]    : string $stable_id
  Example    : my $member = $ma->fetch_by_source_stable_id(
                   "Uniprot/SWISSPROT", "O93279");
  Example    : my $member = $ma->fetch_by_source_stable_id(
                   undef, "O93279");
  Description: Fetches the member corresponding to this $stable_id.
               Although two members from different sources might
               have the same stable_id, this never happens in a normal
               compara DB. You can set the first argument to undef
               like in the second example.
  Returntype : Bio::EnsEMBL::Compara::Member object
  Exceptions : throws if $stable_id is undef
  Caller     : 

=cut

sub fetch_by_source_stable_id {
  my ($self,$source_name, $stable_id) = @_;

  unless(defined $stable_id) {
    throw("fetch_by_source_stable_id must have an stable_id");
  }

  #construct a constraint like 't1.table1_id = 1'
  my $constraint = "";
  $constraint = "m.source_name = '$source_name' AND " if ($source_name);
  $constraint .= "m.stable_id = '$stable_id'";

  #return first element of _generic_fetch list
  my ($obj) = @{$self->_generic_fetch($constraint)};
  return $obj;
}

sub fetch_all_by_source_stable_ids {
  my ($self,$source_name, $stable_ids) = @_;
  return [] if (!$stable_ids or !@$stable_ids);

  #construct a constraint like 't1.table1_id = 1'
  my $constraint = "";
  $constraint = "m.source_name = '$source_name' AND " if ($source_name);
  $constraint .= "m.stable_id IN ('".join("','", @$stable_ids). "')";

  #return first element of _generic_fetch list
  my $obj = $self->_generic_fetch($constraint);
  return $obj;
}

=head2 fetch_all

  Arg        : None
  Example    : my $members = $ma->fetch_all;
  Description: Fetch all the members in the db
  Returntype : listref of Bio::EnsEMBL::Compara::Member objects
  Exceptions : 
  Caller     : 

=cut

sub fetch_all {
  my $self = shift;

  return $self->_generic_fetch();
}


=head2 fetch_by_source

  DEPRECATED: use fetch_all_by_source instead

=cut

sub fetch_by_source {
  my ($self, @args) = @_;
  return $self->fetch_all_by_source(@args);
}

=head2 fetch_all_by_source

  Arg [1]    : string $source_name
  Example    : my $members = $ma->fetch_all_by_source(
                   "Uniprot/SWISSPROT");
  Description: Fetches the member corresponding to a source_name.
  Returntype : listref of Bio::EnsEMBL::Compara::Member objects
  Exceptions : throws if $source_name is undef
  Caller     : 

=cut

sub fetch_all_by_source {
  my ($self,$source_name) = @_;

  throw("source_name arg is required\n")
    unless ($source_name);

  my $constraint = "m.source_name = '$source_name'";

  return $self->_generic_fetch($constraint);
}


=head2 fetch_by_source_taxon

  DEPRECATED: use fetch_all_by_source_taxon instead

=cut

sub fetch_by_source_taxon {
  my ($self, @args) = @_;
  return $self->fetch_all_by_source_taxon(@args);
}

=head2 fetch_all_by_source_taxon

  Arg [1]    : string $source_name
  Arg [2]    : int $taxon_id
  Example    : my $members = $ma->fetch_all_by_source_taxon(
                   "Uniprot/SWISSPROT", 9606);
  Description: Fetches the member corresponding to a source_name and a taxon_id.
  Returntype : listref of Bio::EnsEMBL::Compara::Member objects
  Exceptions : throws if $source_name or $taxon_id is undef
  Caller     : 

=cut

sub fetch_all_by_source_taxon {
  my ($self,$source_name,$taxon_id) = @_;

  throw("source_name and taxon_id args are required") 
    unless($source_name && $taxon_id);

  my $constraint = "m.source_name = '$source_name' and m.taxon_id = $taxon_id";

  return $self->_generic_fetch($constraint);
}

=head2 fetch_all_by_source_genome_db_id

  Arg [1]    : string $source_name
  Arg [2]    : int $genome_db_id
  Example    : my $members = $ma->fetch_all_by_source_genome_db_id(
                   "Uniprot/SWISSPROT", 90);
  Description: Fetches the member corresponding to a source_name and a genome_db_id.
  Returntype : listref of Bio::EnsEMBL::Compara::Member objects
  Exceptions : throws if $source_name or $genome_db_id is undef
  Caller     : 

=cut

sub fetch_all_by_source_genome_db_id {
  my ($self,$source_name,$genome_db_id) = @_;

  throw("source_name and genome_db_id args are required") 
    unless($source_name && $genome_db_id);

  my $constraint = "m.source_name = '$source_name' and m.genome_db_id = $genome_db_id";

  return $self->_generic_fetch($constraint);
}


sub _fetch_all_by_source_taxon_chr_name_start_end_strand_limit {
  my ($self,$source_name,$taxon_id,$chr_name,$chr_start,$chr_end,$chr_strand,$limit) = @_;

  $self->throw("source_name and taxon_id args are required") 
    unless($source_name && $taxon_id && $chr_name && $chr_start && $chr_end && $chr_strand && $limit);

  my $constraint = "m.source_name = '$source_name' and m.taxon_id = $taxon_id 
                    and m.chr_name = \"$chr_name\" 
                    and m.chr_start >= $chr_start and m.chr_end <= $chr_end 
                    and m.chr_strand = $chr_strand limit $limit";

  return $self->_generic_fetch($constraint);
}

sub _fetch_all_by_source_taxon_chr_name_start_end_strand {
  my ($self,$source_name,$taxon_id,$chr_name,$chr_start,$chr_end,$chr_strand) = @_;

  $self->throw("source_name and taxon_id args are required") 
    unless($source_name && $taxon_id && $chr_name && $chr_start && $chr_end && $chr_strand);

  my $constraint = "m.source_name = '$source_name' and m.taxon_id = $taxon_id 
                    and m.chr_name = \"$chr_name\" 
                    and m.chr_start >= $chr_start and m.chr_end <= $chr_end 
                    and m.chr_strand = $chr_strand";

  return $self->_generic_fetch($constraint);
}

=head2 get_source_taxon_count

  Arg [1]    : string $source_name
  Arg [2]    : int $taxon_id
  Example    : my $sp_gene_count = $memberDBA->get_source_taxon_count('ENSEMBLGENE',$taxon_id);
  Description: 
  Returntype : int $sp_gene_count is the number of members for this source_name and taxon_id
  Exceptions : 
  Caller     : 

=cut

sub get_source_taxon_count {
  my ($self,$source_name,$taxon_id) = @_;

  throw("source_name and taxon_id args are required") 
    unless($source_name && $taxon_id);

  my $sth = $self->prepare
    ("SELECT COUNT(*) FROM member WHERE source_name=? AND taxon_id=?");
  $sth->execute($source_name, $taxon_id);
  my ($count) = $sth->fetchrow_array();
  $sth->finish;

  return $count;
}


=head2 fetch_by_relation

  DEPRECATED: use fetch_all_by_relation instead

=cut

sub fetch_by_relation {
  my ($self, @args) = @_;
  return $self->fetch_all_by_relation(@args);
}

=head2 fetch_all_by_relation

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_by_relation {
  my ($self, $relation) = @_;

  my $join;
  my $constraint;

  throw() 
    unless (defined $relation && ref $relation);
  
  if ($relation->isa('Bio::EnsEMBL::Compara::Family')) {
    my $family_id = $relation->dbID;
    $constraint = "fm.family_id = $family_id";
    my $extra_columns = [qw(fm.family_id
                            fm.member_id
                            fm.cigar_line)];
    $join = [[['family_member', 'fm'], 'm.member_id = fm.member_id', $extra_columns]];
  }
  elsif ($relation->isa('Bio::EnsEMBL::Compara::Domain')) {
    my $domain_id = $relation->dbID;
    $constraint = "dm.domain_id = $domain_id";
    my $extra_columns = [qw(dm.domain_id
                            dm.member_id
                            dm.member_start
                            dm.member_end)];
    $join = [[['domain_member', 'dm'], 'm.member_id = dm.member_id', $extra_columns]];
  }
  elsif ($relation->isa('Bio::EnsEMBL::Compara::Homology')) {
    my $homology_id = $relation->dbID;
    $constraint .= "hm.homology_id = $homology_id";
    my $extra_columns = [qw(hm.homology_id
                            hm.member_id
                            hm.peptide_member_id
                            hm.peptide_align_feature_id
                            hm.cigar_line
                            hm.cigar_start
                            hm.cigar_end
                            hm.perc_cov
                            hm.perc_id
                            hm.perc_pos)];
    $join = [[['homology_member', 'hm'], 'm.member_id = hm.member_id', $extra_columns]];
  }
  else {
    throw();
  }

  return $self->_generic_fetch($constraint, $join);
}


=head2 fetch_by_relation_source

  DEPRECATED: use fetch_all_by_relation_source instead

=cut

sub fetch_by_relation_source {
  my ($self, @args) = @_;
  return $self->fetch_all_by_relation_source(@args);
}

=head2 fetch_all_by_relation_source

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_by_relation_source {
  my ($self, $relation, $source_name) = @_;

  throw() 
    unless (defined $relation && ref $relation);
  
  throw("source_name arg is required\n")
    unless ($source_name);

  my $join;
  my $constraint = "m.source_name = '$source_name'";

  if ($relation->isa('Bio::EnsEMBL::Compara::Family')) {
    my $family_id = $relation->dbID;
    $constraint .= " AND fm.family_id = $family_id";
    my $extra_columns = [qw(fm.family_id
                            fm.member_id
                            fm.cigar_line)];
    $join = [[['family_member', 'fm'], 'm.member_id = fm.member_id', $extra_columns]];
  }
  elsif ($relation->isa('Bio::EnsEMBL::Compara::Domain')) {
    my $domain_id = $relation->dbID;
    $constraint .= " AND dm.domain_id = $domain_id";
    my $extra_columns = [qw(dm.domain_id
                            dm.member_id
                            dm.member_start
                            dm.member_end)];
    $join = [[['domain_member', 'dm'], 'm.member_id = dm.member_id', $extra_columns]];
  }
  elsif ($relation->isa('Bio::EnsEMBL::Compara::Homology')) {
    my $homology_id = $relation->dbID;
    $constraint .= " AND hm.homology_id = $homology_id";
    my $extra_columns = [qw(hm.homology_id
                            hm.member_id
                            hm.peptide_member_id
                            hm.peptide_align_feature_id
                            hm.cigar_line
                            hm.cigar_start
                            hm.cigar_end
                            hm.perc_cov
                            hm.perc_id
                            hm.perc_pos)];
    $join = [[['homology_member', 'hm'], 'm.member_id = hm.member_id', $extra_columns]];
  }
  else {
    throw();
  }
  return $self->_generic_fetch($constraint, $join);
}


=head2 fetch_by_relation_source_taxon

  DEPRECATED: use fetch_all_by_relation_source_taxon instead

=cut

sub fetch_by_relation_source_taxon {
  my ($self, @args) = @_;
  return $self->fetch_all_by_relation_source_taxon(@args);
}

=head2 fetch_all_by_relation_source_taxon

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub fetch_all_by_relation_source_taxon {
  my ($self, $relation, $source_name, $taxon_id) = @_;

  throw()
    unless (defined $relation && ref $relation);
  
  throw("source_name and taxon_id args are required") 
    unless($source_name && $taxon_id);

  my $join;
  my $constraint = "m.source_name = '$source_name' AND m.taxon_id = $taxon_id";

  if ($relation->isa('Bio::EnsEMBL::Compara::Family')) {
    my $family_id = $relation->dbID;
    $constraint .= " AND fm.family_id = $family_id";
    my $extra_columns = [qw(fm.family_id
                         fm.member_id
                         fm.cigar_line)];
    $join = [[['family_member', 'fm'], 'm.member_id = fm.member_id', $extra_columns]];
  }
  elsif ($relation->isa('Bio::EnsEMBL::Compara::Domain')) {
    my $domain_id = $relation->dbID;
    $constraint .= " AND dm.domain_id = $domain_id";
    my $extra_columns = [qw(dm.domain_id
                         dm.member_id
                         dm.member_start
                         dm.member_end)];
    $join = [[['domain_member', 'dm'], 'm.member_id = dm.member_id', $extra_columns]];
  }
#  elsif ($relation->isa('Bio::EnsEMBL::Compara::Homology')) {
#  }
  else {
    throw();
  }
  return $self->_generic_fetch($constraint, $join);
}


=head2 fetch_by_subset_id

  DEPRECATED: use fetch_all_by_subset_id instead

=cut

sub fetch_by_subset_id {
  my ($self, @args) = @_;
  return $self->fetch_all_by_subset_id(@args);
}


=head2 fetch_all_by_subset_id

  Arg [1]    : int subset_id
  Example    : @members = @{$memberAdaptor->fetch_all_by_subset_id($subset_id)};
  Description: given a subset_id, does a join to the subset_member table
               to return a list of Member objects in this subset
  Returntype : list by reference of Compara::Member objects
  Exceptions :
  Caller     : general

=cut

sub fetch_all_by_subset_id {
  my ($self, $subset_id) = @_;

  throw() unless (defined $subset_id);

  my $constraint = "sm.subset_id = '$subset_id'";

  my $join = [[['subset_member', 'sm'], 'm.member_id = sm.member_id']];

  return $self->_generic_fetch($constraint, $join);
}


=head2 fetch_gene_for_peptide_member_id

  Arg [1]    : int member_id of a peptide member
  Example    : $geneMember = $memberAdaptor->fetch_gene_for_peptide_member_id($peptide_member_id);
  Description: given a member_id of a peptide member,
               does a join to a copy of member table to extract a member for its gene
  Returntype : Bio::EnsEMBL::Compara::Member object
  Exceptions :
  Caller     : general

=cut

sub fetch_gene_for_peptide_member_id {
  my ($self, $peptide_member_id) = @_;

  throw() unless (defined $peptide_member_id);

  my $constraint = "pepm.member_id = '$peptide_member_id'";

  my $join = [[['member', 'pepm'], 'm.member_id = pepm.gene_member_id']];

  my $obj = undef;
  eval {
    ($obj) = @{$self->_generic_fetch($constraint, $join)};
  };
  return $obj;
}


=head2 fetch_all_peptides_for_gene_member_id

  DEPRECATED: use fetch_all_peptides_for_gene_member_id instead

=cut

sub fetch_peptides_for_gene_member_id {
  my ($self, @args) = @_;
  return $self->fetch_all_peptides_for_gene_member_id(@args);
}

=head2 fetch_all_peptides_for_gene_member_id

  Arg [1]    : int member_id of a gene member
  Example    : @pepMembers = @{$memberAdaptor->fetch_all_peptides_for_gene_member_id($gene_member_id)};
  Description: given a member_id of a gene member,
               fetches all peptide members for this gene
  Returntype : array ref of Bio::EnsEMBL::Compara::Member objects
  Exceptions :
  Caller     : general

=cut

sub fetch_all_peptides_for_gene_member_id {
  my ($self, $gene_member_id) = @_;

  throw() unless (defined $gene_member_id);

  my $constraint = "m.gene_member_id = '$gene_member_id'";

  my $peplist = undef;
  eval {
    $peplist = $self->_generic_fetch($constraint);
  };
  return $peplist;
}


=head2 fetch_canonical_peptide_member_for_gene_member_id

  Arg [1]    : int member_id of a gene member
  Example    : $pepMembers = $memberAdaptor->fetch_peptides_for_gene_member_id($gene_member_id);
  Description: given a member_id of a gene member,
               fetches all peptide members for this gene
  Returntype : Bio::EnsEMBL::Compara::Member object
  Exceptions :
  Caller     : general

=cut

sub fetch_canonical_peptide_member_for_gene_member_id {
  my ($self, $gene_member_id) = @_;

  throw() unless (defined $gene_member_id);

  my $constraint = "m.gene_member_id = '$gene_member_id'";
  my $join = [[['subset_member', 'sm'], 'sm.member_id = m.member_id']];

  #fixed fetch_canonical_peptide_member_for_gene_member_id so that it
  #returns the same canonical peptide used in the
  #peptide_align_feature.  There are some cases where a gene will have
  #multiple transcripts but with the same translation sequence (and
  #hence the same sequence length). The information comes from the
  #ensembl core database and is specified by the canonical_transcript
  #relationship

  my $obj = undef;
  eval {
    ($obj) = @{$self->_generic_fetch($constraint, $join)};
  };
  $self->_final_clause("");
  return $obj;
}


=head2 fetch_canonical_transcript_member_for_gene_member_id

  Arg [1]    : int member_id of a gene member
  Example    : $pepMembers = $memberAdaptor->fetch_peptides_for_gene_member_id($gene_member_id);
  Description: given a member_id of a gene member,
               fetches all peptide members for this gene
  Returntype : Bio::EnsEMBL::Compara::Member object
  Exceptions :
  Caller     : general

=cut

sub fetch_canonical_transcript_member_for_gene_member_id {
  my ($self, $gene_member_id) = @_;

  throw() unless (defined $gene_member_id);

  my $constraint = "m.gene_member_id = '$gene_member_id'";

  my $obj = undef;
  eval {
    ($obj) = @{$self->_generic_fetch($constraint)};
  };
  $self->_final_clause("");
  return $obj;
}


#
# INTERNAL METHODS
#
###################

=head2 _generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) string $logic_name
               the logic_name of the analysis of the features to obtain
  Example    : $fts = $a->_generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Performs a database fetch and returns feature objects in
               contig coordinates.
  Returntype : listref of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : BaseFeatureAdaptor, ProxyDnaAlignFeatureAdaptor::_generic_fetch

=cut
  
sub _generic_fetch {
  my ($self, $constraint, $join) = @_;

  my @tables = @{$self->tables};
  my $columns = join(', ', @{$self->columns()});
  
  if ($join) {
    foreach my $single_join (@{$join}) {
      my ($tablename, $condition, $extra_columns) = @{$single_join};
      if ($tablename && $condition) {
        push @tables, $tablename;
        
        if($constraint) {
          $constraint .= " AND $condition";
        } else {
          $constraint = " $condition";
        }
      } 
      if ($extra_columns) {
        $columns .= ", " . join(', ', @{$extra_columns});
      }
    }
  }
      
  #construct a nice table string like 'table1 t1, table2 t2'
  my $tablenames = join(', ', map({ join(' ', @$_) } @tables));

  my $sql = "SELECT $columns FROM $tablenames";

  my $default_where = $self->_default_where_clause;
  my $final_clause = $self->_final_clause;

  #append a where clause if it was defined
  if($constraint) { 
    $sql .= " WHERE $constraint ";
    if($default_where) {
      $sql .= " AND $default_where ";
    }
  } elsif($default_where) {
    $sql .= " WHERE $default_where ";
  }

  #append additional clauses which may have been defined
  $sql .= " $final_clause" if($final_clause);

  #print("$sql\n");
  my $sth = $self->prepare($sql);
  $sth->execute;

#  print STDERR $sql,"\n";
  return $self->_objs_from_sth($sth);
}

sub tables {
  return [['member', 'm']];
}

sub columns {
  return ['m.member_id',
          'm.source_name',
          'm.stable_id',
          'm.version',
          'm.taxon_id',
          'm.genome_db_id',
          'm.description',
          'm.chr_name',
          'm.chr_start',
          'm.chr_end',
          'm.chr_strand',
          'm.sequence_id',
          'm.gene_member_id',
          'm.display_label'
          ];
}

sub create_instance_from_rowhash {
	my ($self, $rowhash) = @_;
	
	return Bio::EnsEMBL::Compara::Member->new_fast({
		_dbID           => $rowhash->{member_id},
		_stable_id      => $rowhash->{stable_id},
		_version        => $rowhash->{version},
		_taxon_id       => $rowhash->{taxon_id},
		_genome_db_id   => $rowhash->{genome_db_id},
		_description    => $rowhash->{description},
		_chr_name       => $rowhash->{chr_name},
		_chr_start      => $rowhash->{chr_start},
		_chr_end        => $rowhash->{chr_end},
		_chr_strand     => $rowhash->{chr_strand},
		_sequence_id    => $rowhash->{sequence_id} || 0,
		_source_name    => $rowhash->{source_name},
		_display_label  => $rowhash->{display_label},
		_gene_member_id => $rowhash->{gene_member_id},
		_adaptor        => $self
	});
}

sub init_instance_from_rowhash {
  my $self = shift;
  my $member = shift;
  my $rowhash = shift;

  $member->member_id($rowhash->{'member_id'});
  $member->stable_id($rowhash->{'stable_id'});
  $member->version($rowhash->{'version'});
  $member->taxon_id($rowhash->{'taxon_id'});
  $member->genome_db_id($rowhash->{'genome_db_id'});
  $member->description($rowhash->{'description'});
  $member->chr_name($rowhash->{'chr_name'});
  $member->chr_start($rowhash->{'chr_start'});
  $member->chr_end($rowhash->{'chr_end'});
  $member->chr_strand($rowhash->{'chr_strand'});
  $member->sequence_id($rowhash->{'sequence_id'});
  $member->gene_member_id($rowhash->{'gene_member_id'});
  $member->source_name($rowhash->{'source_name'});
  $member->display_label($rowhash->{'display_label'});
  $member->adaptor($self);

  return $member;
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my @members = ();

  while(my $rowhash = $sth->fetchrow_hashref) {
    my ($member,$attribute);
    $member = $self->create_instance_from_rowhash($rowhash);
    
    my @_columns = @{$self->columns};
    if (scalar keys %{$rowhash} > scalar @_columns) {
      $attribute = new Bio::EnsEMBL::Compara::Attribute;
      $attribute->member_id($rowhash->{'member_id'});
      foreach my $autoload_method (keys %$rowhash) {
        next if (grep /$autoload_method/,  @_columns);
        $attribute->$autoload_method($rowhash->{$autoload_method});
      }
    }
    if (defined $attribute) {
      push @members, [$member, $attribute];
    } else {
      push @members, $member;
    } 
  }
  $sth->finish;
  return \@members
}

sub _default_where_clause {
  my $self = shift;
  return '';
}

sub _final_clause {
  my $self = shift;

  $self->{'_final_clause'} = shift if(@_);
  return $self->{'_final_clause'};
}

sub _fetch_sequence_by_id {
  my ($self, $sequence_id) = @_;
  return $self->db->get_SequenceAdaptor->fetch_by_dbID($sequence_id);
}

sub _fetch_sequence_exon_bounded_by_member_id {
  my ($self, $member_id) = @_;
  return $self->db->get_SequenceAdaptor->fetch_sequence_exon_bounded_by_member_id($member_id);
}

sub _fetch_sequence_cds_by_member_id {
  my ($self, $member_id) = @_;
  return $self->db->get_SequenceAdaptor->fetch_sequence_cds_by_member_id($member_id);
}


sub create_AlignedMember_from_member_attribute {
  my $self = shift;
  my $member_attribute = shift;

  my ($gene_member, $attribute) = @{$member_attribute};
  my $member = $self->fetch_by_dbID($attribute->peptide_member_id);

  bless $member, "Bio::EnsEMBL::Compara::AlignedMember";
  $member->cigar_line($attribute->cigar_line);
  $member->cigar_start($attribute->cigar_start);
  $member->cigar_end($attribute->cigar_end);
  $member->adaptor(undef);

  return $member;
}



#
# STORE METHODS
#
################

=head2 store

  Arg [1]    : 
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : 

=cut

sub store {
  my ($self,$member) = @_;

  unless($member->isa('Bio::EnsEMBL::Compara::Member')) {
    throw(
      "member arg must be a [Bio::EnsEMBL::Compara::Member]"
    . "not a $member");
  }

  my $sth = $self->prepare("INSERT ignore INTO member (stable_id,version, source_name,
                              taxon_id, genome_db_id, description,
                              chr_name, chr_start, chr_end, chr_strand,display_label)
                            VALUES (?,?,?,?,?,?,?,?,?,?,?)");

  my $insertCount = $sth->execute($member->stable_id,
                  $member->version,
                  $member->source_name,
                  $member->taxon_id,
                  $member->genome_db_id,
                  $member->description,
                  $member->chr_name,
                  $member->chr_start,
                  $member->chr_end,
                  $member->chr_strand,
                  $member->display_label);
  if($insertCount>0) {
    #sucessful insert
    $member->dbID( $sth->{'mysql_insertid'} );
    $sth->finish;
  } else {
    $sth->finish;
    #UNIQUE(source_name,stable_id) prevented insert since member was already inserted
    #so get member_id with select
    my $sth2 = $self->prepare("SELECT member_id, sequence_id FROM member WHERE source_name=? and stable_id=?");
    $sth2->execute($member->source_name, $member->stable_id);
    my($id, $sequence_id) = $sth2->fetchrow_array();
    warn("MemberAdaptor: insert failed, but member_id select failed too") unless($id);
    $member->dbID($id);
    $member->sequence_id($sequence_id) if ($sequence_id);
    $sth2->finish;
  }

  $member->adaptor($self);

  # insert in sequence table to generate new
  # sequence_id to insert into member table;
  if(defined($member->sequence) and $member->sequence_id == 0) {
    $member->sequence_id($self->db->get_SequenceAdaptor->store($member->sequence,1)); # Last parameter induces a check for redundancy

    my $sth3 = $self->prepare("UPDATE member SET sequence_id=? WHERE member_id=?");
    $sth3->execute($member->sequence_id, $member->dbID);
    $sth3->finish;
  }

  return $member->dbID;
}

sub store_reused {
  my ($self,$member) = @_;

  unless($member->isa('Bio::EnsEMBL::Compara::Member')) {
    throw(
      "member arg must be a [Bio::EnsEMBL::Compara::Member]"
    . "not a $member");
  }

  my $sth = $self->prepare("INSERT ignore INTO member (member_id, stable_id, version, source_name,
                              taxon_id, genome_db_id, description,
                              chr_name, chr_start, chr_end, chr_strand,display_label)
                            VALUES (?,?,?,?,?,?,?,?,?,?,?,?)");

  my $insertCount = $sth->execute(
                  $member->member_id,
                  $member->stable_id,
                  $member->version,
                  $member->source_name,
                  $member->taxon_id,
                  $member->genome_db_id,
                  $member->description,
                  $member->chr_name,
                  $member->chr_start,
                  $member->chr_end,
                  $member->chr_strand,
                  $member->display_label);
  if($insertCount>0) {
    #sucessful insert
    $member->dbID( $sth->{'mysql_insertid'} );
    $sth->finish;
  } else {
    $sth->finish;
    #UNIQUE(source_name,stable_id) prevented insert since member was already inserted
    #so get member_id with select
    my $sth2 = $self->prepare("SELECT member_id, sequence_id FROM member WHERE source_name=? and stable_id=?");
    $sth2->execute($member->source_name, $member->stable_id);
    my($id, $sequence_id) = $sth2->fetchrow_array();
    warn("MemberAdaptor: insert failed, but member_id select failed too") unless($id);
    $member->dbID($id);
    $member->sequence_id($sequence_id) if ($sequence_id);
    $sth2->finish;
  }

  $member->adaptor($self);

  # insert in sequence table to generate new
  # sequence_id to insert into member table;
  if(defined($member->sequence) and $member->sequence_id == 0) {
    $member->sequence_id($self->db->get_SequenceAdaptor->store($member->sequence,1)); # Last parameter induces a check for redundancy

    my $sth3 = $self->prepare("UPDATE member SET sequence_id=? WHERE member_id=?");
    $sth3->execute($member->sequence_id, $member->dbID);
    $sth3->finish;
  }

  return $member->dbID;
}


sub update_sequence {
  my ($self, $member) = @_;

  return 0 unless($member);
  unless($member->dbID) {
    throw("MemberAdapter::update_sequence member must have valid dbID\n");
  }
  unless(defined($member->sequence)) {
    warning("MemberAdapter::update_sequence with undefined sequence\n");
  }

  if($member->sequence_id) {
    my $sth = $self->prepare("UPDATE sequence SET sequence = ?, length=? WHERE sequence_id = ?");
    $sth->execute($member->sequence, $member->seq_length, $member->sequence_id);
    $sth->finish;
  } else {
    $member->sequence_id($self->db->get_SequenceAdaptor->store($member->sequence,1)); # Last parameter induces a check for redundancy

    my $sth3 = $self->prepare("UPDATE member SET sequence_id=? WHERE member_id=?");
    $sth3->execute($member->sequence_id, $member->dbID);
    $sth3->finish;
  }
  return 1;
}

sub store_gene_peptide_link {
  my ($self, $gene_member_id, $peptide_member_id) = @_;

  eval {
    my $sth = $self->prepare("UPDATE member SET gene_member_id=? where member_id=?");
    $sth->execute($gene_member_id, $peptide_member_id);
    $sth->finish;
  };
}

# DEPRECATED METHODS
####################

sub fetch_longest_peptide_member_for_gene_member_id {
  my $self = shift;

  throw("Method deprecated. You can now use the fetch_canonical_peptide_member_for_gene_member_id method\n");
}

1;

