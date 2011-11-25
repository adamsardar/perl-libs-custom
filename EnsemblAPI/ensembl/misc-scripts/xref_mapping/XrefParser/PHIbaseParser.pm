package XrefParser::PHIbaseParser;

use strict;
use POSIX qw(strftime);
use File::Basename;

use XML::LibXML;

use base qw( XrefParser::BaseParser );

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print "\nUsage: PhibaseParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  print STDERR "ARGV[0]: $ARGV[0]\n";

  run($ARGV[0]);

}

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id    = shift;
  my $species_id   = shift;
  #my $phi_xml_file = shift;
  my $files        = shift;
  my $verbose      = shift;

  my $phi_xml_file = @{$files}[0];
  
  print STDERR "PhiBase file to parse, $phi_xml_file\n";

  my %phi_mapping;
  my %taxIds;

  my $term = undef;
  my $desc = undef;
  
  my $phi_parser = XML::LibXML->new();
  my $phi_doc    = $phi_parser->parse_file($phi_xml_file);

  my $concepts_count = 0;
  
  foreach my $concept ($phi_doc->findnodes('ondex:ondexdata/ondexdataseq/concepts/concept')) {

      #print STDERR "parsing concept...\n";

      my ($pid_node) = $concept->findnodes('./pid');
      my $pid = $pid_node->to_literal;

      #print STDERR "pid, $pid\n";

    my $concept_accessions_aref = $concept->findnodes('./coaccessions/concept_accession');
    my $uniprot_acc = undef;
    foreach my $concept_accession (@$concept_accessions_aref) {
	
	# get the one associated with UniProt

	my ($elementOf) = $concept_accession->findnodes('./elementOf');

	# print STDERR "elementOf: " . $elementOf->to_literal . "\n";

	if ($elementOf->to_literal =~ /UPROT/) {
	    my ($accession) = $concept_accession->findnodes('./accession');
	    $uniprot_acc = $accession->to_literal;

	    #print STDERR "UniProt: $uniprot_acc\n";
	    
	}
    }

    if (!defined $uniprot_acc) {
	# print STDERR "phi id, $pid, no uniprot mapping found in xml file!\n";
    }
    else {
	if (!defined $phi_mapping{$uniprot_acc}) {
	    $phi_mapping{$uniprot_acc} = 
		{
		    -phi_ids => [$pid],
		    -tax_id  => undef,
		};
	}
	else {
	    my $phi_href = $phi_mapping{$uniprot_acc};
	    my $aref = $phi_href->{-phi_ids};
	    push (@$aref, $pid);
	}
    }

      # Get the TaxId

      my $concept_gds_aref = $concept->findnodes('./cogds/concept_gds');
      my $taxId = undef;
      
      foreach my $concept_gds (@$concept_gds_aref) {
	  
	  # get the one associated with Taxid
	  
	  my ($attrname) = $concept_gds->findnodes('./attrname');
	  
	  if ($attrname->to_literal =~ /TAXID/) {
	      my ($value) = $concept_gds->findnodes('./value');
	      $taxId = $value->to_literal;
	      $taxId =~ s/\D//g;
	      
	      #print STDERR "taxId: $taxId\n";

	      if (! defined $taxIds{$taxId}) {
		  $taxIds{$taxId} = 1;
	      }

	      if (defined $uniprot_acc) {
		  my $phi_href = $phi_mapping{$uniprot_acc};
		  $phi_href->{-tax_id} = $taxId;
	      }
	  }
      }
      
      #print STDERR "parsing concept done\n\n";
      
      $concepts_count++;
  }
  
  my @phis = keys (%phi_mapping);

  print STDERR "Parsed $concepts_count concepts\n";
  print STDERR "Found " . @phis . " with UniProt mapping!\n";

  print STDERR "Found " . keys (%taxIds) . " different taxIds\n";
  
  #get the "main" PHIbase source id.
  $source_id = $self->get_source_id_for_source_name("PHIbase");


  #get the mapping that are already there so that we don't get lots of duplicates.
  # stored in the global hash xref_dependent_mapped.
  $self->get_dependent_mappings($source_id);

  #if(!defined($species_id)){
  #  $species_id = $self->get_species_id_for_filename($phi_xml_file);
  #}

  my $swiss_miss=0;
  my (%swiss) = %{$self->get_valid_codes("uniprot/", $species_id)};

  print STDERR "got " . keys (%swiss) . " Uniprot entries\n";
   
  print STDERR "species_id, source_id: $species_id, $source_id\n";

  # Don't check only the species_id, but all taxIds specified in xref_config.ini

  my %species2tax = $self->species_id2taxonomy();
  my @tax_ids = @{$species2tax{$species_id}};

  print STDERR "tax_ids from xref_config.ini file: " . join (', ', @tax_ids) . "\n";

  my $added = 0;

  foreach my $uniprot_acc (keys (%phi_mapping)) {
      my $phis_href = $phi_mapping{$uniprot_acc};
      my $taxId = $phis_href->{-tax_id};
      if (grep {$_ eq $taxId} @tax_ids) {

	  print STDERR "Adding xrefs for UniProt, $uniprot_acc\n";

	  # Get the master_xref_id
	  # and the linkage
	  
	  my $master_xref_id = $swiss{$uniprot_acc};

	  if (!defined $master_xref_id) {
	      print STDERR "failed to get the master_xref_if for UniProt, $uniprot_acc!\n";
	      # one reason it happens is that the UniProt identifier is attached to a different tax node in Phibase that it is in UniProt
	  }
	  else {
	      print STDERR "master_xref_id, $master_xref_id\n";
	      
	      my $linkage = undef;
	      
	      my $phis_aref = $phis_href->{-phi_ids};
	      foreach my $phibase_id (@$phis_aref) {
		  print STDERR "Adding xrefs for phibase id, $phibase_id\n";
		  $self->add_to_xrefs($master_xref_id,$phibase_id,'',$phibase_id,'',$linkage,$source_id,$species_id);
		  $added++;
	      }
	  }

      }
  }

  print STDERR "Added $added PHIbase xrefs\n";

  return 0;

}

1;
