package Bio::EnsEMBL::Compara::DBSQL::ConstrainedElementAdaptor;

use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Compara::ConstrainedElement;
use Bio::EnsEMBL::Compara::DnaFrag;
use Bio::EnsEMBL::Utils::Exception;
use Data::Dumper;

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  $self->{_use_autoincrement} = 0;
  return $self;
}

=head2 store

  Arg  1     : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object 
  Arg  2     : listref of Bio::EnsEMBL::Compara::ConstrainedElement ($constrained_element) objects 
               The things you want to store
  Example    : none
  Description: It stores the given ConstrainedElements in the database.
  Returntype : none
  Exceptions : throw if Arg-1 is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object
	       throw if Arg-2 is not a Bio::EnsEMBL::Compara::ConstrainedElement	
  Caller     : called by the Bio::EnsEMBL::Compara::Production::GenomicAlignBlock::Gerp module 

=cut

sub store {
    my ( $self, $mlss_obj, $constrained_elements ) = @_;
    if (defined($mlss_obj)) {
	throw("$mlss_obj is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object")	
	  unless ($mlss_obj->isa("Bio::EnsEMBL::Compara::MethodLinkSpeciesSet"));
    } else {
	throw("undefined Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object");
    }
    
    my $mlssid = $mlss_obj->dbID;
    
    #
    #Find unique constrained_element_id by using a temporary table with an auto_increment column
    #
    my $ce_id_sql = "INSERT INTO constrained_element_production (constrained_element_id) VALUES (\"NULL\")";
    my $ce_id_sth = $self->prepare($ce_id_sql);

    my $constrained_element_sql = qq{INSERT INTO constrained_element (
		constrained_element_id,
		dnafrag_id,
		dnafrag_start, 
		dnafrag_end,
		dnafrag_strand,
		score,
		method_link_species_set_id,
		p_value,
		taxonomic_level
	) VALUES (?,?,?,?,?,?,?,?,?)};
    
    my $constrained_element_sth = $self->prepare($constrained_element_sql) or die;
    
    foreach my $constrained_element_group (@$constrained_elements) {
	$ce_id_sth->execute();
	my $constrained_element_id = $ce_id_sth->{'mysql_insertid'};
	if ($constrained_element_id < $mlssid * 10000000000 || 
	    $constrained_element_id > ($mlssid+1) * 10000000000) {
	    $constrained_element_id = $mlssid * 10000000000 + $constrained_element_id;
	}

	foreach my $constrained_element (@{$constrained_element_group}) {
	    throw("$constrained_element is not a Bio::EnsEMBL::Compara::ConstrainedElement object")
	      unless ($constrained_element->isa("Bio::EnsEMBL::Compara::ConstrainedElement"));
	    $constrained_element_sth->execute(
					      $constrained_element_id,
					      $constrained_element->reference_dnafrag_id,
					      $constrained_element->start,
					      $constrained_element->end,
					      $constrained_element->strand,
					      $constrained_element->score,
					      $mlssid,
					      ($constrained_element->p_value or undef),
					      ($constrained_element->taxonomic_level or undef)
					     );
	}
    }
}

=head2 delete_by_MethodLinkSpeciesSet

  Arg  1     : method_link_species_set object $mlss
  Example    : $constrained_element_adaptor->delete_by_MethodLinkSpeciesSet($mlss);
  Description: It removes constrained elements with the specified method_link_species_set_id from the database
  Returntype : none
  Exceptions : throw if passed parameter is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object 
  Caller     : general

=cut

sub delete_by_MethodLinkSpeciesSet {
  my ($self, $mlss_obj) = @_;

  if (defined($mlss_obj)) {
    throw("$mlss_obj is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object")      
      unless ($mlss_obj->isa("Bio::EnsEMBL::Compara::MethodLinkSpeciesSet"));
  } else {
    throw("undefined Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object");
  }

  my $cons_ele_sql =
        qq{DELETE FROM constrained_element WHERE method_link_species_set_id = ?};
  
# Delete constrtained element entries by mlss_id
  my $sth = $self->prepare($cons_ele_sql);
  $sth->execute($mlss_obj->dbID);
  $sth->finish;
}


=head2 delete_by_dbID

  Arg  1     : int $constrained_element_id
  Example    : $constrained_element_adaptor->delete_by_dbID(123);
  Description: It removes constrained elements with the specified ID
  Returntype : none
  Exceptions : none
  Caller     : general

=cut

sub delete_by_dbID {
  my ($self, $constrained_element_id) = @_;

  if (!defined($constrained_element_id)) {
    throw("undefined Constrained Element ID");
  }

  my $cons_ele_sql =
        qq{DELETE FROM constrained_element WHERE constrained_element_id = ?};
  
# Delete constrtained element entries by mlss_id
  my $sth = $self->prepare($cons_ele_sql);
  $sth->execute($constrained_element_id);
  $sth->finish;
}


=head2 fetch_all_by_MethodLinkSpeciesSet_Slice

  Arg  1     : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet $mlss_obj
  Arg  2     : Bio::EnsEMBL::Slice $slice_obj
  Example    : my $listref_of_constrained_elements = $constrained_element_adaptor->
  		fetch_all_by_MethodLinkSpeciesSet_Slice($mlss_obj, $slice_obj);
  Description: Retrieve the corresponding
               Bio::EnsEMBL::Compara::ConstrainedElement object listref
  Returntype : listref of Bio::EnsEMBL::Compara::ConstrainedElement objects
  Exceptions : throw if Arg-1 is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object
	       throw if Arg-2 is not a Bio::EnsEMBL::Slice object
  Caller     : object::methodname

=cut

sub fetch_all_by_MethodLinkSpeciesSet_Slice {
	my ($self, $mlss_obj, $slice_obj) = @_;
	if (defined($mlss_obj)) {
		throw("$mlss_obj is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object")	
		unless ($mlss_obj->isa("Bio::EnsEMBL::Compara::MethodLinkSpeciesSet"));
	} else {
		throw("undefined Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object");
	}
	if (defined($slice_obj)) {
		throw("$slice_obj is not a Bio::EnsEMBL::Slice object")	
		unless ($slice_obj->isa("Bio::EnsEMBL::Slice"));
	} else {
		throw("undefined Bio::EnsEMBL::Slice object");
	}

	my $dnafrag_adp = $self->db->get_DnaFragAdaptor;
	my $dnafrag = $dnafrag_adp->fetch_by_Slice($slice_obj);
	my $sql = qq{
		WHERE
		method_link_species_set_id = ?
		AND
		dnafrag_id = ? 
	};
	my (@constrained_elements, $lower_bound);

	if(defined($slice_obj->start) && defined($slice_obj->end) && 
		($slice_obj->start <= $slice_obj->end)) {
			my $max_alignment_length = $mlss_obj->max_alignment_length;
			$lower_bound = $slice_obj->start - $max_alignment_length;
			$sql .= qq{
				AND
				dnafrag_end >= ?
				AND
				dnafrag_start <= ?
				AND
				dnafrag_start >= ?
			};
	}
	
	$self->_fetch_all_ConstrainedElements($sql, \@constrained_elements,
		$mlss_obj->dbID, $dnafrag->dbID, $slice_obj->start, $slice_obj->end, $lower_bound, $slice_obj);
	return \@constrained_elements;
}

=head2 fetch_all_by_MethodLinkSpeciesSet_Dnafrag

  Arg  1     : Bio::EnsEMBL::Compara::MethodLinkSpeciesSet mlss_obj
  Arg  2     : Bio::EnsEMBL::Compara::DnaFrag dnafrag_obj
  Example    : my $listref_of_constrained_elements = $constrained_element_adaptor->
  		fetch_all_by_MethodLinkSpeciesSet_Dnafrag($mlss_obj, $dnafrag_obj);
  Description: Retrieve the corresponding
               Bio::EnsEMBL::Compara::ConstrainedElement object listref
  Returntype : listref of Bio::EnsEMBL::Compara::ConstrainedElement objects
  Exceptions : throw if Arg-1 is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet mlss_obj object
	       throw if Arg-2 is not a Bio::EnsEMBL::Compara::DnaFrag object
  Caller     : object::methodname

=cut

sub fetch_all_by_MethodLinkSpeciesSet_Dnafrag {
	my ($self, $mlss_obj, $dnafrag_obj, $dnafrag_start, $dnafrag_end) = @_;
	if (defined($mlss_obj)) {
		throw("$mlss_obj is not a Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object")	
		unless ($mlss_obj->isa("Bio::EnsEMBL::Compara::MethodLinkSpeciesSet"));
	} else {
		throw("undefined Bio::EnsEMBL::Compara::MethodLinkSpeciesSet object");
	} 
	if(defined($dnafrag_obj)) {
		throw("$dnafrag_obj is not a Bio::EnsEMBL::Compara::DnaFrag object")
		unless ($dnafrag_obj->isa("Bio::EnsEMBL::Compara::DnaFrag"));
	} else {
		throw("undefined Bio::EnsEMBL::Compara::DnaFrag object");
	}

	my (@constrained_elements, $lower_bound);
	my $sql = qq{
		WHERE
		method_link_species_set_id = ?
		AND
		dnafrag_id = ?
	};
 
	if (defined($dnafrag_start) && defined($dnafrag_end) && ($dnafrag_start <= $dnafrag_end)) {
		my $max_alignment_length = $mlss_obj->max_alignment_length;
		$lower_bound = $dnafrag_start - $max_alignment_length;
	} else {
		$dnafrag_start = 1;
		$dnafrag_end = $dnafrag_obj->length;
		$lower_bound = $dnafrag_start;
	}
	$sql .= qq{
		AND
		dnafrag_end >= ?
		AND
		dnafrag_start <= ?
		AND
		dnafrag_start >= ?
	};
	$self->_fetch_all_ConstrainedElements($sql, \@constrained_elements, 
			$mlss_obj->dbID, $dnafrag_obj->dbID, $dnafrag_start, $dnafrag_end, $lower_bound);
	return \@constrained_elements;
}

sub _fetch_all_ConstrainedElements {#used when getting constrained elements by slice or dnafrag
	my ($self) = shift;
	my ($sql, $constrained_elements, $mlss_id, $dnafrag_id, $start, $end, $lower_bound, $slice) = @_;
	$sql = qq{
       		SELECT
       		constrained_element_id,
       		dnafrag_start,
       		dnafrag_end,
                dnafrag_strand,
      		score,
      		p_value,
       		taxonomic_level
       		FROM
       		constrained_element} . $sql;

	my $sth = $self->prepare($sql);
	$sth->execute($mlss_id, $dnafrag_id, $start, $end, $lower_bound);
	my ($dbID, $ce_start, $ce_end, $ce_strand, $score, $p_value, $tax_level);
	$sth->bind_columns(\$dbID, \$ce_start, \$ce_end, \$ce_strand, \$score, \$p_value, \$tax_level);
	while ($sth->fetch()) {
		my $constrained_element = Bio::EnsEMBL::Compara::ConstrainedElement->new_fast (
			{
				'adaptor' => $self,
				'dbID' => $dbID,
				'slice' => $slice,
				'start' =>  ($ce_start - $start + 1), 
				'end' => ($ce_end - $start + 1),
			        'strand' => $ce_strand,
				'method_link_species_set_id' => $mlss_id,
				'score' => $score,
				'p_value' => $p_value,
				'taxonomic_level' => $tax_level,
				'reference_dnafrag_id' => $dnafrag_id,
			}
		);
		push(@$constrained_elements, $constrained_element);
	}
}	

=head2 fetch_all_by_dbID

  Arg  1     : listref of constrained_element_ids
  Example    : my $listref_of_constrained_elements = $constrained_element_adaptor->fetch_all_by_dbID($list_ref_of_constrained_element_ids);
  Description: Retrieve the corresponding constrained_elements from a given list of constrained_element_ids 
  Returntype : listref of Bio::EnsEMBL::Compara::ConstrainedElement constrained_elements 
  Exceptions : throw if Arg-1 is not a listref
  Caller     : object::methodname

=cut

sub fetch_all_by_dbID {
	my ($self, $constrained_element_ids) = @_;
	if(defined($constrained_element_ids)) {
		throw("Arg-1 needs to be a listref of dbIDs") unless (
			ref($constrained_element_ids) eq "ARRAY");
	}
	my @constrained_elements;
	my $sql = qq{
		WHERE
		ce.constrained_element_id = ?
	};
	$self->_fetch_all_ConstrainedElements_by_dbID($sql, \@constrained_elements, $constrained_element_ids);
	return \@constrained_elements;
}

=head2 fetch_by_dbID

  Arg  1     : int constrained_element_id
  Example    : my $constrained_element = $constrained_element_adaptor->
               fetch_by_dbID($constrained_element_id);
  Description: Retrieve the corresponding constrained_element.
  Returntype : Bio::EnsEMBL::Compara::ConstrainedElement object
  Exceptions : -none-
  Caller     : object::methodname

=cut

sub fetch_by_dbID {
  my ($self, $constrained_element_id) = @_;
  return ($self->fetch_all_by_dbID([$constrained_element_id]))->[0];
}

sub _fetch_all_ConstrainedElements_by_dbID {#used when getting constrained elements by constrained_element_id
	my ($self) = shift;
	my ($sql, $constrained_elements, $dbIDs) = @_;
        
	$sql = qq{
       		SELECT
       		ce.constrained_element_id,
       		ce.dnafrag_id,
       		ce.dnafrag_start,
       		ce.dnafrag_end,
       		ce.method_link_species_set_id,
      		ce.score,
      		ce.p_value,
       		ce.taxonomic_level,
		gdb.name,
		df.name
       		FROM
       		constrained_element ce
		INNER JOIN 
		dnafrag df 
		ON
		df.dnafrag_id = ce.dnafrag_id
		INNER JOIN 
		genome_db gdb
		ON 
		gdb.genome_db_id = df.genome_db_id} . $sql;

	my $sth = $self->prepare($sql);
	foreach my $constrained_element_id (@{ $dbIDs }) {
		my (%general_attributes, @alignment_segments);
		$sth->execute( $constrained_element_id );
		my ($dbID, $dnafrag_id, $ce_start, $ce_end, $mlssid, $score, $p_value, $tax_level, $species_name, $dnafrag_name);
	 	$sth->bind_columns(\$dbID, \$dnafrag_id, \$ce_start, \$ce_end, \$mlssid, 
					\$score, \$p_value, \$tax_level, \$species_name, \$dnafrag_name);	
		while ($sth->fetch()) {
			$general_attributes{dbID} = $dbID;
			$general_attributes{mlssid} = $mlssid;
			$general_attributes{score} = $score;
			$general_attributes{p_value} = $p_value;
			$general_attributes{taxonomic_level} = $tax_level;
			push(@alignment_segments, [ $dnafrag_id, $ce_start, $ce_end, $species_name, $dnafrag_name ]);
		}
		my $constrained_element = Bio::EnsEMBL::Compara::ConstrainedElement->new_fast (
			{
				'adaptor' => $self,
				'dbID' => $general_attributes{dbID},
				'alignment_segments' => \@alignment_segments,
				'method_link_species_set_id' => $general_attributes{mlssid},
				'score' => $general_attributes{score},
				'p_value' => $general_attributes{p_value},
				'taxonomic_level' => $general_attributes{taxonomic_level},
			}
		);
		push(@$constrained_elements, $constrained_element);
	}
}

1;
