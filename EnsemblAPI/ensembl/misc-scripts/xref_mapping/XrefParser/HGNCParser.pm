package XrefParser::HGNCParser;

use strict;
use File::Basename;

use base qw( XrefParser::BaseParser );

my $xref_sth ;
my $dep_sth;
my $syn_sth;

# --------------------------------------------------------------------------------
# Parse command line and run if being run directly

if (!defined(caller())) {

  if (scalar(@ARGV) != 1) {
    print STDERR "\nUsage: HGNCParser.pm file <source_id> <species_id>\n\n";
    exit(1);
  }

  run(@ARGV);
}

sub run {

  my $self = shift if (defined(caller(1)));

  my $source_id = shift;
  my $species_id = shift;
  my $files_ref  = shift;
  my $rel_file   = shift;
  my $verbose = shift;

  my $file = @{$files_ref}[0];

  if(!defined($source_id)){
    $source_id = XrefParser::BaseParser->get_source_id_for_filename($file);
  }
  if(!defined($species_id)){
    $species_id = XrefParser::BaseParser->get_species_id_for_filename($file);
  }


  my $hgnc_refseq_manual = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","refseq_manual");
  if(!defined($hgnc_refseq_manual)){
    die  "Could not get source id for HGNC with priority description of refseq_manual\n";
  }
  my $hgnc_refseq_mapped = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","refseq_mapped");
  if(!defined($hgnc_refseq_mapped)){
    die  "Could not get source id for HGNC with priority description of refseq_mapped\n";
  }
  my $hgnc_swissprot_manual = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","swissprot_manual");
  if(!defined($hgnc_swissprot_manual)){
    die  "Could not get source id for HGNC with priority description of swissprot_manual\n";
  }

  my $hgnc_entrezgene_manual  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","entrezgene_manual");
  if(!defined($hgnc_entrezgene_manual)){
    die  "Could not get source id for HGNC with priority description of entrezgene_manual\n";
  }
  my $hgnc_entrezgene_mapped  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","entrezgene_mapped");
  if(!defined($hgnc_entrezgene_mapped)){
    die  "Could not get source id for HGNC with priority description of entrezgene_mapped\n";
  }

  my $hgnc_ensembl_mapped  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","ensembl_manual");
  if(!defined($hgnc_ensembl_mapped)){
    die  "Could not get source id for HGNC with priority description of ensembl_manual\n";
  }

  my $hgnc_desc_only  = XrefParser::BaseParser->get_source_id_for_source_name("HGNC","desc_only");
  if(!defined($hgnc_desc_only)){
    die  "Could not get source id for HGNC with priority description of desc_only\n";
  }  

  my $hgnc_lrg  = XrefParser::BaseParser->get_source_id_for_source_name("LRG_HGNC_notransfer");
  if(!defined($hgnc_lrg)){
    die  "Could not get source id for LRG_HGNC_notransfer\n";
  }

  my (%swissprot)  =  %{XrefParser::BaseParser->get_valid_codes("Uniprot/SWISSPROT",$species_id)};
  my (%refseq) =  %{XrefParser::BaseParser->get_valid_codes("refseq",$species_id)};
  my @list;
  push @list, "refseq_peptide";
  push @list, "refseq_dna";
  my (%entrezgene) = %{XrefParser::BaseParser->get_valid_xrefs_for_dependencies("EntrezGene",@list)};

  my $refseq_count = 0;
  my $swissprot_count = 0;
  my $entrezgene_count = 0;
  my $ensembl_count = 0;
  my $mismatch = 0;

  my $hugo_io = $self->get_filehandle($file);

  if ( !defined $hugo_io ) {
    print "ERROR: Can't open HGNC file $file\n";
    return 1;
  }

  $_ = $hugo_io->getline();

  while ( $_ = $hugo_io->getline() ) {

    chomp;

    # 0 HGNC ID	           # primary accession
    # 1 Approved Symbol    # label
    # 2 Approved Name      # description
    # 3 Previous Symbols   # synonyms
    # 4 Aliases            # aliases
    # 5 entrezgene ID   manually curated
    # 6 RefSeq ID       manually curated
    # 7 entrezgene ID   mapped
    # 8 RefSeq ID       mapped
    # 9 Ensembl ID     manual
    #10 Swissprot ID  manual
    #11 <strong>.....</strong>  look for lrg specification

    my @array = split(/\t/,$_);

    # Use the RefSeq if available as this is manually curated
    # If no RefSeq, use the Swissprot instead

    my $seen = 0;
    
    if($array[11] and $array[11] =~ /http:\/\/www.lrg-sequence.org\/LRG\/(LRG_\d+)\'/){
#      print $array[0]." linked to LRG $1\n";
      my $lrg_stable_id = $1;
      XrefParser::BaseParser->add_to_direct_xrefs($lrg_stable_id, 'gene', $array[0], '', $array[1], $array[2], "", $hgnc_lrg, $species_id);
      
      $self->add_synonyms_for_hgnc($hgnc_lrg, $array[0], $species_id, $array[3], $array[4]);
    }
    
    if ($array[9]){              # Ensembl direct xref
      $seen = 1;
      $ensembl_count++;
      XrefParser::BaseParser->add_to_direct_xrefs($array[9],'gene', $array[0], '', $array[1], $array[2], "", $hgnc_ensembl_mapped, $species_id);
      
      $self->add_synonyms_for_hgnc($hgnc_ensembl_mapped, $array[0], $species_id, $array[3], $array[4]);
      
    }
    if ($array[6]) {             # RefSeq
      if(defined($refseq{$array[6]})){
	$seen = 1;
	$refseq_count++;
	XrefParser::BaseParser->add_to_xrefs($refseq{$array[6]}, $array[0], '', $array[1], $array[2], "", $hgnc_refseq_manual, $species_id);
	$self->add_synonyms_for_hgnc($hgnc_refseq_manual, $array[0], $species_id, $array[3], $array[4]);
      }	
    }
    if ($array[10]) {             # Swissprot
      if(defined($swissprot{$array[10]})){
	$seen = 1;
	$swissprot_count++;
	XrefParser::BaseParser->add_to_xrefs($swissprot{$array[10]}, $array[0], '', $array[1], $array[2], "", $hgnc_swissprot_manual, $species_id);
	$self->add_synonyms_for_hgnc($hgnc_swissprot_manual, $array[0], $species_id, $array[3], $array[4]);
	
      }
    }	
    if ($array[8]) {             # RefSeq
      if(defined($refseq{$array[8]})){
	$seen = 1;
	$refseq_count++;
	XrefParser::BaseParser->add_to_xrefs($refseq{$array[8]}, $array[0], '', $array[1], $array[2], "", $hgnc_refseq_mapped, $species_id);
	
	$self->add_synonyms_for_hgnc($hgnc_refseq_mapped, $array[0], $species_id, $array[3], $array[4]);
      }
    }
    
    if(defined($array[5])){
      if(defined($entrezgene{$array[5]})){
	$seen = 1;
	XrefParser::BaseParser->add_to_xrefs($entrezgene{$array[5]}, $array[0], '', 
					     $array[1], $array[2], "", $hgnc_entrezgene_manual, $species_id);
	$entrezgene_count++;
	$self->add_synonyms_for_hgnc($hgnc_entrezgene_manual, $array[0], $species_id, $array[3], $array[4]);
	
      }
    }
    
    if(defined($array[7])){
      if(defined($entrezgene{$array[7]})){
	$seen = 1;
	XrefParser::BaseParser->add_to_xrefs($entrezgene{$array[7]}, $array[0], '', 
					     $array[1], $array[2], "", $hgnc_entrezgene_mapped, $species_id);
	$entrezgene_count++;
	$self->add_synonyms_for_hgnc($hgnc_entrezgene_mapped, $array[0], $species_id, $array[3], $array[4]);
      }    
    }
    if(!$seen){ # Store to keep descriptions etc
      $self->add_xref($array[0], "", $array[1], $array[2], $hgnc_desc_only, $species_id, "MISC");      
      $self->add_synonyms_for_hgnc($hgnc_desc_only, $array[0], $species_id, $array[3], $array[4]);
      $mismatch++;
    }
  }


  $hugo_io->close();
  
  print "Loaded a total of " . ($refseq_count + $entrezgene_count + $swissprot_count) . " HGNC xrefs, $refseq_count from RefSeq curated mappings and $entrezgene_count from EntrezGene mappings $swissprot_count from Swissprot and $ensembl_count from ensembl_mapping\n" if($verbose);
  
  print "$mismatch xrefs could not be associated via RefSeq, EntrezGene or ensembl\n" if($verbose);
  
  return 0; # successful
  
}


sub add_synonyms_for_hgnc{
  my $self = shift;
  my $type = shift;
  my $name = shift;
  my $species_id = shift;
  my $dead_name = shift;
  my $alias = shift;

    
    if (defined($dead_name)) {     # dead name, add to synonym
      my @array2 = split(',\s*', $dead_name);
      foreach my $arr (@array2){
	XrefParser::BaseParser->add_to_syn($name, $type, $arr, $species_id);
      }
    }
  
  if (defined($alias)) {     # alias, add to synonym
    my @array2 = split(',\s*', $alias);
    foreach my $arr (@array2){
      XrefParser::BaseParser->add_to_syn($name, $type, $arr, $species_id);
    }
  }
}


sub rename_url_file{
  return "hugo.txt";
}

1;
    

