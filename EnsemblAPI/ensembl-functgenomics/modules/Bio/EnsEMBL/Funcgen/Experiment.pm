
=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Experiment
  
=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Experiment;

my $array = Bio::EnsEMBL::Funcgen::Experiment->new(
						   -ADAPTOR             => $self,
						   -NAME                => $name,
					       -GROUP               => $group,
						   -DATE                => $date,
						   -PRIMARY_DESIGN_TYPE => $p_design_type,
						   -DESCRIPTION         => $description,
                                                   );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $exp_adaptor = $db_adaptor->get_ExperimentAdaptor();
my $exp = $exp_adaptor->fetch_by_name($exp_name)

=head1 DESCRIPTION

An Experiment object represents an experiment instance . The data
are stored in the experiment, egroup, target, design_type and 
experimental_variable tables.

=cut


################################################################################

package Bio::EnsEMBL::Funcgen::Experiment;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Storable;


use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Funcgen::Storable);



=head2 new

  Arg [-NAME]: string - the name of this experiment
  Arg [-GROUP]: string - the name of the experimental group
  Arg [-GROUP_ID]: int - the dbID of the experimental group
  Arg [-DATE]: string - the date of the experiment (format?)
  Arg [-PRIMARY_DESIGN_TYPE]: string - MGED term for the primary design of teh experiment e.g. binding_site_identification
  Arg [-DESCRIPTION]: string - of the experiment

  Example    : my $array = Bio::EnsEMBL::Funcgen::Experiment->new(
								  -NAME                => $name,
								  -GROUP               => $group,
								  -DATE                => $date,
								  -PRIMARY_DESIGN_TYPE => $p_design_type,
								  -DESCRIPTION         => $description,
                                                 		 );
  Description: Creates a new Bio::EnsEMBL::Funcgen::Experiment object.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : None ? Should throw if mandatory params not set
  Caller     : General
  Status     : Medium Risk

=cut

#experimental_variables?
#design_type
#target(s)?

sub new {
	my $caller = shift;

	my $class = ref($caller) || $caller;

	my $self = $class->SUPER::new(@_);

	my ($name, $group_id, $group, $date, $p_dtype, $desc, $xml_id, $xml)
		= rearrange( ['NAME', 'GROUP_ID', 'GROUP', 'DATE', 'PRIMARY_DESIGN_TYPE', 'DESCRIPTION', 'MAGE_XML', 'MAGE_XML_ID'], @_ );
	
	$self->name($name)          if defined $name;
	$self->group_id($group_id)  if defined $group_id;
	$self->group($group)        if defined $group;
	$self->date($date)          if defined $date;
	$self->primary_design_type($p_dtype)    if defined $p_dtype;
	$self->description($desc)   if defined $desc;
	$self->mage_xml_id($xml_id) if defined $xml_id;
	$self->mage_xml($xml)       if defined $xml;


	#Need to add mandatory params check here!!
	#name, group or group_id
	

	return $self;
}


### GENERIC ACCESSOR METHODS ###

=head2 name

  Arg [1]: string - the name of this experiment
  Example: $exp->name('Experiment-1');
  Description: Getter/Setter for the experiment name
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name{
	my ($self) = shift;	

	$self->{'name'} = shift if(@_);

	return $self->{'name'};
}

=head2 group_id

  Arg [1]: int - the group dbID corresponding to this experiment
  Example: $exp->group_db_id('1');
  Description: Getter/Setter for the group_db_id
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut



sub group_id{
	my ($self) = shift;	

	$self->{'group_id'} = shift if(@_);

	return $self->{'group_id'};
}

=head2 mage_xml

  Arg [1]: string(optional) - MAGE XML
  Example: my $xml = $exp->mage_xml();
  Description: Getter/Setter for the mage_xml attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut



sub mage_xml{
  my ($self) = shift;	

  $self->{'mage_xml'} = shift if(@_);
  
  #use fetch_attrs?
  if(! exists $self->{'mage_xml'} && $self->mage_xml_id()){
	$self->{'mage_xml'} = $self->adaptor->fetch_mage_xml_by_Experiment($self);
  }

  return (exists $self->{'mage_xml'}) ? $self->{'mage_xml'} : undef;
}

=head2 mage_xml_id

  Arg [1]: int (optional) - mage_xml_id
  Example: $exp->group_db_id('1');
  Description: Getter/Setter for the mage_xml attribute
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut



sub mage_xml_id{
  my $self = shift;	

  $self->{'mage_xml_id'} = shift if @_;
  
  return $self->{'mage_xml_id'};
}





=head2 group

  Arg [1]: optional - Bio::EnsEMBL::Funcgen::Group
  Example: my $exp_group_name = $exp->group->name();
  Description: Getter/Setter for the group
  Returntype : Bio::EnsEMBL::Funcgen::Group
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut


sub group{
  my ($self) = shift;	

  if(@_){
    $self->{'group'} = shift;
  }elsif(! exists $self->{'group'}){
    #this should be group_adaptor->fetch_by_dbID($self->group_id());
    throw("Not yet implemented group adaptor");
  }

  return $self->{'group'};
}

=head2 date

  Arg [1]: optional - date, format yyyy-mm-dd
  Example: $exp->date('2006-06-09');
  Description: Getter/Setter for the date
  Returntype : date string
  Exceptions : None ? should throw/warn if format not correct
  Caller     : General
  Status     : Medium

=cut

sub date{
  my $self = shift;
  
  if(@_){
    #Need to validate format here
    $self->{'date'} = shift;
  }

  return $self->{'date'};
}

=head2 description

  Arg [1]: string - the experiment description
  Example: $exp->description("Human chromosome X TFBS identification");
  Description: Getter/Setter for the experiment description
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub description{
  my $self = shift;
  $self->{'description'} = shift if(@_);
  return $self->{'description'};
}

=head2 primary_design_type

  Arg [1]: string - MGED term for primary design type
  Example: $exp->primary_design_type('binding_site_identification');
  Description: Getter/Setter for the primary design type
  Returntype : string
  Exceptions : None ? should throw if not MGED term
  Caller     : General
  Status     : At risk

=cut

sub primary_design_type{
  my ($self) = shift;
	
  if(@_){
    #warn "Need to validate design_types against MGED";
    $self->{'primary_design_type'} = shift;
  }
  return $self->{'primary_design_type'};
}



#These convenience methods are to provide a registry for the experimental chips of the experiment

=head2 get_ExperimentalChips

  Example: my $exp_chips = @{$exp->get_ExperimentalChips()}
  Description: Retrieves all ExperiemntalChips
  Returntype : Listref of ExperimentalChips
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub get_ExperimentalChips{
  my ($self) = shift;
	
  #should this also store echips?

  #Need to retrieve all from DB if not defined, then check whether already present and add and store if not
  #what about immediate access to dbID
  #should we separate and have add_ExperimentalChip?

  

  if(! exists $self->{'experimental_chips'}){
     $self->{'experimental_chips'} = {};

     #need to warn about DBAdaptor here?
  
    foreach my $echip(@{$self->adaptor->db->get_ExperimentalChipAdaptor->fetch_all_by_experiment_dbID($self->dbID())}){
      $self->{'experimental_chips'}->{$echip->unique_id()} = $echip;
    }
  }

  #is this returning a list or a listref?
  return [values %{$self->{'experimental_chips'}}];
}

=head2 add_ExperimentalChip

  Example     : $exp_chip = $exp->add_ExperimentalChip($exp_chip)
  Description : Adds and stores an ExperiemntalChip for this Experiment
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions  : Throws is not passed a valid stored Bio::EnsENBML::Funcgen::ExperimentalChip
  Caller      : General
  Status      : At risk

=cut

sub add_ExperimentalChip{
  my ($self, $echip) = @_;
 

 throw("Must pass a valid stored Bio::EnsEMBL::Funcgen::ExperimentalChip object") 
    if(! $echip->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip") || ! $echip->dbID());
 
  if(! exists $self->{'experimental_chips'}){
    $self->get_ExperimentalChips();
    $self->{'experimental_chips'}->{$echip->unique_id()} = $echip;
    #do this here without checking to avoid probelm of retrieving first stored chip
  }elsif(exists  $self->{'experimental_chips'}->{$echip->unique_id()}){
    warn("You cannot add the same ExperimentalChip(".$echip->unique_id().")to an Experiment more than once, check your code");
  }else{
    $self->{'experimental_chips'}->{$echip->unique_id()} = $echip;
  }
  
  return;
}

=head2 get_ExperimentalChip_by_unique_id

  Example: $exp_chip = $exp->add_ExperimentalChip($exp_chip)
  Description: Adds and stores an ExperiemntalChip for this Experiment
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Exceptions : Throws if no uid supplied
  Caller     : General
  Status     : At risk

=cut

sub get_ExperimentalChip_by_unique_id{
  my ($self, $uid) = @_;
  
  my ($echip);

  throw("Must supply a ExperimentalChip unque_id") if(! defined $uid);
  
  $self->{'experimental_chips'} || $self->get_ExperimentalChips();

  if(exists $self->{'experimental_chips'}->{$uid}){
    $echip = $self->{'experimental_chips'}->{$uid};
  }
  #should we warn here if not exists?

  return $echip;
}


=head2 get_ExperimentalChip_unique_ids

  Example:     foreach my $uid(@{$self->experiment->get_ExperimentalChip_unique_ids()}){ ... }
  Description: retrieves all ExperimentalChip unique_ids
  Returntype : ListRef
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub get_ExperimentalChip_unique_ids{
  my $self = shift;
  
  $self->{'experimental_chips'} || $self->get_ExperimentalChips();

  return [keys %{ $self->{'experimental_chips'}}];
}



#should we add a methods to return just the 


#methods?
#lazy load design_types and exp_variables
#target?  Is this a one to one?



1;

