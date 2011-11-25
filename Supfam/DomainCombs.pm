package Supfam::DomainCombs;
require Exporter;

=head1 NAME

Supfam::DomainCombs.pm

=head1 SYNOPSIS

Package for investigating SUPERFAMILY domain combinations.

=head1 AUTHOR

Matt Oates (Matt.Oates@bristol.ac.uk)

=head1 COPYRIGHT

Copyright 2010 Gough Group, University of Bristol.

=head1 SEE ALSO

Supfam::SQLFunc.pm - Where all the SQL related basic functions are kept.

=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
                     getGenomeUniqDomCombs
                     calcDomPairs
                     getGenomeDomCombs
					removeSharedDomCombs
					splitIntoUPSpeciesCombs
					splitIntoNCBISpeciesCombs
					splitIntoNCBISpeciesCombIDs
					splitIntoNCBISpeciesCombsAbundance
					splitIntoNCBISpeciesFamsAbundance
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use Data::Dumper;
use Term::ProgressBar;
use Math::Combinatorics;
use Supfam::SQLFunc;
use Supfam::Utils;

=pod
=item * getGenomeUniqDomCombs($genome_id)
Returns a hashref of all the unique superfamily combinations and their frequency
for a given genome id.
=cut
sub getGenomeUniqDomCombs {
my ($genome, $dbh) = @_;
$dbh = Supfam::SQLFunc::dbConnect() unless defined $dbh;
my $close_dbh = (@_ < 2)?1:0;

	my $uniq_combs = {};
	getGenomeDomCombs($genome,$uniq_combs,$dbh);
	removeSharedDomCombs($genome,$uniq_combs,$dbh);
	dbDisconnect($dbh) if $close_dbh;

	return $uniq_combs;
}

=pod
=item * getUniqDomCombsExclTaxon($genome_id, $taxon_id)
Returns a hashref of all the unique superfamily combinations for a given genome id.
Uniqueness is calculated in respect to excluding comparissons with all genomes 
below the specified NCBI taxonomy id.
=cut
sub getUniqDomCombsExclTaxon($$) {
}

=pod
=item * getUniqCombsExclGenomes($genome_id, $genomes = [])
Returns a hashref of all the unique superfamily combinations for a given genome id.
Uniqueness is calculated in respect to excluding comparissons with all genomes specified
in the $genomes arrayref.
=cut
sub getUniqDomCombsExclGenomes($\@) {
my ($genome, $to_exclude) = @_;
ref $to_exclude eq "ARRAY" or die "Expected an ARRAY ref $!";
}

=pod
= item * calcDomPairs(pointer to a hash of $->{Domain1}{Domain2}=Frequency of domain pair,pointer to list (array)
of domain archtectures). Returns a pointer to the same hash as before, updated so as to include all of the 
domain pairs present in input list.
 
=back

=cut

sub calcDomPairs($$){
	
	my ($DomPairsHash, $ArchitecturesArray) = @_;
	
	 foreach my $comb (@$ArchitecturesArray) {
                my @combinations = combine(2,grep(!/_gap_/,split(/,/,$comb)));
                foreach my $pair (@combinations) {
                        @_ = sort {$a <=> $b} @$pair;
                        $DomPairsHash->{shift()}{shift()}++;
                }
        }
        return $DomPairsHash;
}

=pod
=item * getGenomeDomCombs($genome, $PointerToCombsHash, $dbh - optional)
Given a genome id and a pointer to a hash (structure: Dom1 => Dom2 => count), this function will
connect to superfamily (or another DB if a database handle is given) and extract the domain combinations
found in that genome. These are then broken down into pairs (excluding _gap_ entries) and pushed onto the
hash pointer. Returned is the hash pointer.
=cut

sub getGenomeDomCombs {
my ($genome, $combs, $dbh) = @_;
$dbh = dbConnect() unless defined $dbh;
my $close_dbh = (@_ < 3)?1:0;

	   my $query = $dbh->prepare("SELECT comb FROM len_comb WHERE genome = ?");
        $query->execute($genome) or return undef;
    
        while (my ($comb) = $query->fetchrow_array()) {
                my @combinations = combine(2,grep(!/_gap_/,split(/,/,$comb)));
                foreach my $pair (@combinations) {
                        @_ = sort {$a <=> $b} @$pair;
                        $combs->{shift()}{shift()}++;
                }
        }
	dbDisconnect($dbh) if $close_dbh;
	return $combs;
}

=pod
=item * removeSharedDomCombs($genome, $PointerToCombsHash, $dbh - optional)
Given a genome id and a pointer to a hash (structure: Dom1 => Dom2 => count), this function will
connect to superfamily (or another DB if a database handle is given) and remove all domain pairs 
found in both the genome of interest and the domain pair hash. Useful for identifying sets of unique
domain pairs.
=cut

sub removeSharedDomCombs {
my ($genome, $combs, $dbh) = @_;
$dbh = dbConnect() unless defined $dbh;
my $close_dbh = (@_ < 3)?1:0;

	my ($nrows) = $dbh->selectrow_array("SELECT count(genome) FROM len_comb WHERE genome != ?", undef, $genome);

	my $query = $dbh->prepare("SELECT comb FROM len_comb WHERE genome != ?");
	$query->execute($genome) or return undef;
	my $pbar = Term::ProgressBar->new({'name' => "Removing shared combs for $genome",
                                           'count' => $nrows,
                                           'remove' => 1,
                                           'ETA' => 'linear',
                                           'fh' => \*STDERR
                                        });
	$pbar->minor(0);

	my $progress = 1;
	my $update = 0;
        while (my ($comb) = $query->fetchrow_array()) {
                my @combinations = combine(2,grep(!/_gap_/,split(/,/,$comb)));
                foreach my $pair (@combinations) {
                        @_ = sort {$a <=> $b} @$pair;
			#Delete this specific dom pair
                        delete $combs->{shift()}{shift()};
                }
                $progress++;
		$update = $pbar->update($progress) if $progress >= $update;
        }
	#Delete hanging empty keys where we removed all pair partners
	while (my ($key,$val) = each(%$combs)) {delete $combs->{$key} unless keys %$val;}
	dbDisconnect($dbh) if $close_dbh;
}

=pod

=item * splitIntoUPSpeciesCombs($genome_id,$SpeciescombsHashRef)
for use with the larger uniprot psuedo-genomes in SUPERFAMILY (e.g. v9 and up). Returns a hashref
 of all the species in the input genome and their domain archtectures (hash=> species => [list of domain
 architectures]).
 
 Works well with uniprot viral ('v9') and whole uniprot ('up') data. Works by looking in the protein table
 for the species name in the 'OS= *******' field.
 
 The two progress bars are there because this can be a LONG job if you are trying to do it for all of Uniprot
 'up')
=cut

sub splitIntoUPSpeciesCombs{
	
	no warnings 'uninitialized';	
	
	my ($genome, $SpeciesCombHashRef, $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	
	my $close_dbh = (@_ < 3)?1:0;
	
	my $query = $dbh->prepare("select count(protein)  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $no_elements = $query->fetchrow_array();
	$query -> finish;
	
    $query = $dbh->prepare("select protein, comment  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $ProteinCommentHash = {};
	
	my $pbar = Term::ProgressBar->new({'name' => "Getting proteins in $genome",
                                           'count' => $no_elements,
                                           'remove' => 1,
                                           'ETA' => 'linear',
                                           'fh' => \*STDERR
                                        });    
        $pbar->minor(0);   
        my $progress = 1;
        my $update = 0;
	
	while (my (@DBresult) = $query->fetchrow_array()) {
	
	$ProteinCommentHash -> {$DBresult[0]}=$DBresult[1];
	#Hash of protein => comment.
	
	$update = $pbar->update($progress) if $progress >= $update;
    $progress++;  
	}
	
	$pbar->update($no_elements);

	   my $nrows = keys(%$ProteinCommentHash);
       $pbar = Term::ProgressBar->new({'name' => "Getting combs for proteins in $genome",
                                           'count' => $nrows,
                                           'remove' => 1,
                                           'ETA' => 'linear',
                                           'fh' => \*STDERR
                                        });    
        $pbar->minor(0);                         
	   
        $progress = 1;
        $update = 0;
	
	while (my ($protein, $comment) = each(%$ProteinCommentHash)){
		
		$query = $dbh->prepare_cached("select comb from comb where protein = ?");
		$query->execute($protein) or return undef;
		
		my $comb = $query->fetchrow_array();

		$comment =~ m/OS=(.*?)[\.]?( [A-Z]{2}=|$)/;
		my $species = $1;
		die "No Match (this script isn't doing what you want it to) - check input for 'OS=': $comment" unless($species);
			
		$SpeciesCombHashRef -> {$species} = [] unless ($SpeciesCombHashRef -> {$species});
		push(@{$SpeciesCombHashRef -> {$species}},$comb);
				
		$query->finish;
		
	    $update = $pbar->update($progress) if $progress >= $update;
		$progress++;   
	}
	
	$pbar->update($nrows);
	
	dbDisconnect($dbh) if $close_dbh;
	return $SpeciesCombHashRef;
}

=pod

=item * splitIntoNCBISpeciesCombs($genome_id,$SpeciescombsHashRef)
for use with the larger NCBI psuedo-genomes in SUPERFAMILY (e.g. vl). Returns a hashref
 of all the species in the input genome and their domain archtectures (hash= species => [list of domain
 architectures]). Known to work on viral NCBI data ('vl').
 
 It works by looking at the comments section of the protein table and extracting te species name
 in [] brackets. As there are sometimes multiple entries in [] brackets, each one is then checked against the NCBI
 names list in Superfamily. Only entries found in that table are included. 
 splitIntoNCBISpeciesCombIDs
 If the species that you're looking for are stored differenty' make your own version of this subroutine and change the regex.
 
=cut

sub splitIntoNCBISpeciesCombs{
	
	no warnings 'uninitialized';	
	
	my ($genome, $SpeciesCombHashRef, , $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	
	my $close_dbh = (@_ < 3)?1:0;
	
	my $NCBINames = {};
	#This is simply a list of all the valid NCBI names in the NCBI treee. This is going to be used as a check later when we pull species names from the comments in the protein table
	
	if(-e './.NCBINames.dat'){
		
		$NCBINames = EasyUnDump('./.NCBINames.dat');
	}else{
		
		my $sth = $dbh->prepare("select name from ncbi_names") or die"Query failed";
		$sth->execute() or return undef;
		
		while(my $name = $sth->fetchrow_array()){
			$NCBINames->{lc($name)}=1;
			#lc (lower case) was used because some entries in ncbi_names were not of the same case as those in the Protein comment table	
		}
		
		EasyDump('./.NCBINames.dat',$NCBINames);
	}
	
	my $query = $dbh->prepare("select count(protein)  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $no_elements = $query->fetchrow_array();
	$query -> finish;
	
    $query = $dbh->prepare("select protein, comment  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $ProteinCommentHash = {};
	
	while (my (@DBresult) = $query->fetchrow_array()) {
	
	$ProteinCommentHash -> {$DBresult[0]}=$DBresult[1];
	#Hash of protein => comment.
	}
	
	   my $nrows = keys(%$ProteinCommentHash);
       my $pbar = Term::ProgressBar->new({'name' => "Getting combs for proteins in $genome",
                                           'count' => $nrows,
                                           'remove' => 1,
                                           'ETA' => 'linear',
                                           'fh' => \*STDERR
                                        });    
        $pbar->minor(0);                         
	   
        my $progress = 1;
        my $update = 0;
	
	while (my ($proteinid, $comment) = each(%$ProteinCommentHash)){
		
		$query = $dbh->prepare_cached("select comb from comb where protein = ?");
		$query->execute($proteinid) or return undef;
		
		my $comb = $query->fetchrow_array();
		
		my @SpeciesTags;
		
		while($comment =~ s/\[(.*?\[.*\].*?|.*?)\]//){
	
			push(@SpeciesTags,$1);
		}
		
		
		my $flag =0; #This is just a quick check to make sure that at least one of the comments corresponds to an NCBI specie		
			
		foreach my $species (@SpeciesTags){
			
			if($NCBINames->{lc($species)}){
				
				$SpeciesCombHashRef -> {$species} = [] unless ($SpeciesCombHashRef -> {$species});
				push(@{$SpeciesCombHashRef -> {$species}},$comb);
					
				$flag = 1;
				}
			}
			
		die "\nNo NCBI species were found in the comment for protein: $proteinid Specie: @SpeciesTags  \n" unless($flag);
			
		$query->finish;
		
	    $update = $pbar->update($progress) if $progress >= $update;
		$progress++;   
	}
	
	$pbar->update($nrows);
	
	dbDisconnect($dbh) if $close_dbh;
	return $SpeciesCombHashRef; # A hash 0f $Hash->{species}[domain combinations]
}

=pod

=item * splitIntoNCBISpeciesCombIDs($genome_id,$SpeciescombsHashRef, $fresh)
for use with the larger NCBI psuedo-genomes in SUPERFAMILY (e.g. vl). Returns a hashref
 of all the species in the input genome and their domain archtectures (hash= species => [list of domain
 combinations IDs]). Known to work on viral NCBI data ('vl').
 
 It works by looking at the comments section of the protein table and extracting te species name
 in [] brackets. As there are sometimes multiple entries in [] brackets, each one is then checked against the NCBI
 names list in Superfamily. Only entries found in that table are included. 
 
 If the species that you're looking for are stored differenty' make your own version of this subroutine and change the regex.
 
=cut

sub splitIntoNCBISpeciesCombIDs{
	
	no warnings 'uninitialized';	
	
	my ($genome, $SpeciesCombHashRef, , $fresh, $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	
	my $close_dbh = (@_ < 4)?1:0;
	
	my $NCBINames = {};
	#This is simply a list of all the valid NCBI names in the NCBI treee. This is going to be used as a check later when we pull species names from the comments in the protein table
	
	if(-e './.NCBINames.dat' && ! $fresh){
		
		$NCBINames = EasyUnDump('./.NCBINames.dat');
	}else{
		
		my $sth = $dbh->prepare("select name from ncbi_names") or die"Query failed";
		$sth->execute() or return undef;
		
		while(my $name = $sth->fetchrow_array()){
			$NCBINames->{lc($name)}=1;
			#lc (lower case) was used because some entries in ncbi_names were not of the same case as those in the Protein comment table	
		}
		
		EasyDump('./.NCBINames.dat',$NCBINames);
	}
	
	my $query = $dbh->prepare("select count(protein)  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $no_elements = $query->fetchrow_array();
	$query -> finish;
	
    $query = $dbh->prepare("select protein, comment  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $ProteinCommentHash = {};
	
	while (my (@DBresult) = $query->fetchrow_array()) {
	
	$ProteinCommentHash -> {$DBresult[0]}=$DBresult[1];
	#Hash of protein => comment.
	}
	
	   my $nrows = keys(%$ProteinCommentHash);
       my $pbar = Term::ProgressBar->new({'name' => "Getting combs for proteins in $genome",
                                           'count' => $nrows,
                                           'remove' => 1,
                                           'ETA' => 'linear',
                                           'fh' => \*STDERR
                                        });    
        $pbar->minor(0);                         
	   
        my $progress = 1;
        my $update = 0;
	
	while (my ($proteinid, $comment) = each(%$ProteinCommentHash)){
		
		$query = $dbh->prepare_cached("select comb_id from comb where protein = ?");
		$query->execute($proteinid) or return undef;
		
		my $combid = $query->fetchrow_array();
		
		my @SpeciesTags;
		
		while($comment =~ s/\[(.*?\[.*\].*?|.*?)\]//){
	
			push(@SpeciesTags,$1);
		}
		
		
		my $flag =0; #This is just a quick check to make sure that at least one of the comments corresponds to an NCBI specie		
			
		foreach my $species (@SpeciesTags){
			
			if($NCBINames->{lc($species)}){
				
				$SpeciesCombHashRef -> {$species} = [] unless ($SpeciesCombHashRef -> {$species});
				push(@{$SpeciesCombHashRef -> {$species}},$combid);
					
				$flag = 1;
				}
			}
			
		die "\nNo NCBI species were found in the comment for protein: $proteinid Specie: @SpeciesTags  \n" unless($flag);
			
		$query->finish;
		
	    $update = $pbar->update($progress) if $progress >= $update;
		$progress++;   
	}
	
	$pbar->update($nrows);
	
	dbDisconnect($dbh) if $close_dbh;
	return $SpeciesCombHashRef; # A hash 0f $Hash->{species}[domain combination IDs]
}

=pod

=item * splitIntoNCBISpeciesFamsAbundance($genome_id,$SpeciesFamsHashRef)
for use with the larger NCBI psuedo-genomes in SUPERFAMILY (e.g. vl). Returns a hashref
 of all the species in the input genome and their families (hash= species => {protein_family_ids = abundance}). Known to work on viral NCBI data ('vl').
 
 It works by looking at the comments section of the protein table and extracting the species name
 in [] brackets. As there are sometimes multiple entries in [] brackets, each one is then checked against the NCBI
 names list in Superfamily. Only entries found in that table are included. 
 
 If the species that you're looking for are stored differenty' make your own version of this subroutine and change the regex.
  
=cut

sub splitIntoNCBISpeciesFamsAbundance{
	
	##This is not my finest piece of perl. But it is faily robust => Adam
	
	my ($genome, $SpeciesFamsHashRef, $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	
	my $close_dbh = (@_ < 3)?1:0;
	
	my $NCBINames = {};
	#This is simply a list of all the valid NCBI names in the NCBI treee. This is going to be used as a check later when we pull species names from the comments in the protein table
	
	my $sth = $dbh->prepare("select name from ncbi_names") or die"Query failed";
	$sth->execute() or return undef;
		
	while(my $name = $sth->fetchrow_array()){
		$NCBINames->{lc($name)}=1;
		#lc (lower case) was used because some entries in ncbi_names were not of the same case as those in the Protein comment table - URGH
	}
	
    my $query = $dbh->prepare("select protein, comment  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $ProteinCommentHash = {};
	
	while (my ($Protein, $Comment) = $query->fetchrow_array()) {
	
	$ProteinCommentHash -> {$Protein}=$Comment;
	#Hash of protein => comment.
	}
	
	$query->finish;
	
	#A hash of proteins and their comments is prepared seperately to parsing the comments field and assigning a species to it. this is due to the horrible way that SF stores the NCBI genomes	
	
	my $FamiliesInDataset = {};	
	
	while (my ($Protein, $comment) = each(%$ProteinCommentHash)){
		
		no warnings 'uninitialized';
		
		$query = $dbh->prepare_cached("SELECT fa FROM ass JOIN family ON ass.auto=family.auto where ass.protein = ?");
		$query->execute($Protein) or return undef;
		
		my $family = $query->fetchrow_array();

		$family = 'Unassigned' if ($family =~ m/^$/);
		
		my @SpeciesTags;
		
		while($comment =~ s/\[(.*?\[.*\].*?|.*?)\]//){
	
			push(@SpeciesTags,$1);
		} #This is all the possible species that this proteins could belong to - we shall compare them against the NCBI table to see which is the 'correct' one (note, allows for more than one, depending on whether the species name is in the NCBInames list we created earlier)
		
		
		my $flag =0; #This is just a quick check to make sure that at least one of the comments corresponds to an NCBI specie		
			
		foreach my $species (@SpeciesTags){
			
			if($NCBINames->{lc($species)}){
				
				$SpeciesFamsHashRef -> {$species} = {} unless ($SpeciesFamsHashRef -> {$species});
				$SpeciesFamsHashRef -> {$species}{$family}++;
				
			    $FamiliesInDataset->{$family}++;
				
				$flag = 1;
				}
		}
			
		die "\n\nNo NCBI species were found in the comment for protein: $Protein Specie: @SpeciesTags  !!!!\n\n" unless($flag);
			
	$query->finish;
	}
		
	dbDisconnect($dbh) if $close_dbh;
	
	return ($FamiliesInDataset); # A hash 0f $Hash->{species}{protein_family}=number of times seen.
}

=item * splitIntoNCBISpeciesCombsAbundance($genome_id,$SpeciesCombsHashRef)
for use with the larger NCBI psuedo-genomes in SUPERFAMILY (e.g. vl). Returns a hashref
 of all the species in the input genome and their families (hash= species => {protein_comb_ids = abundance}). Known to work on viral NCBI data ('vl').
 
 It works by looking at the comments section of the protein table and extracting the species name
 in [] brackets. As there are sometimes multiple entries in [] brackets, each one is then checked against the NCBI
 names list in Superfamily. Only entries found in that table are included. 
 
 If the species that you're looking for are stored differenty' make your own version of this subroutine and change the regex.
  
=cut

sub splitIntoNCBISpeciesCombsAbundance{
	
	##This is not my finest piece of perl. But it is faily robust => Adam
	
	no warnings 'uninitialized';	
	
	my ($genome, $SpeciesCombHashRef, $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	
	my $close_dbh = (@_ < 3)?1:0;
	
	my $NCBINames = {};
	#This is simply a list of all the valid NCBI names in the NCBI treee. This is going to be used as a check later when we pull species names from the comments in the protein table
	
	my $sth = $dbh->prepare("select name from ncbi_names") or die"Query failed";
	$sth->execute() or return undef;
		
	while(my $name = $sth->fetchrow_array()){
		$NCBINames->{lc($name)}=1;
		#lc (lower case) was used because some entries in ncbi_names were not of the same case as those in the Protein comment table - URGH
	}
	
    my $query = $dbh->prepare("select protein, comment  from protein where genome = ?") or die"Query failed";
	$query->execute($genome) or return undef;
	
	my $ProteinCommentHash = {};
	
	while (my ($Protein, $Comment) = $query->fetchrow_array()) {
	
	$ProteinCommentHash -> {$Protein}=$Comment;
	#Hash of protein => comment.
	}

	
	$query->finish;
	
	#A hash of proteins and their comments is prepared seperately to parsing the comments field and assigning a species to it. this is due to the horrible way that SF stores the NCBI genomes	
	
	my $CombsInDataset = {};	
	
	while (my ($Protein, $comment) = each(%$ProteinCommentHash)){
		
		$query = $dbh->prepare_cached("SELECT id FROM comb_index JOIN comb ON comb_index.id=comb.comb_id where comb.protein = ?");
		$query->execute($Protein) or return undef;
		
		my $CombID = $query->fetchrow_array();
		
		my @SpeciesTags;
		
		while($comment =~ s/\[(.*?\[.*\].*?|.*?)\]//){
	
			push(@SpeciesTags,$1);
		} #This is all the possible species that this proteins could belong to - we shall compare them against the NCBI table to see which is the 'correct' one (note, allows for more than one, depending on whether the species name is in the NCBInames list we created earlier)
		
		
		my $flag =0; #This is just a quick check to make sure that at least one of the comments corresponds to an NCBI specie		
			
		foreach my $species (@SpeciesTags){
			
			if($NCBINames->{lc($species)}){
				
				$SpeciesCombHashRef -> {$species} = {} unless ($SpeciesCombHashRef -> {$species});
				$SpeciesCombHashRef -> {$species}{$CombID}++;
				
			    $CombsInDataset->{$CombID}++; #Easy way to count the total number of items in dataset
				
				$flag = 1;
				}
		}
			
		die "\n\nNo NCBI species were found in the comment for protein: $Protein Specie: @SpeciesTags  !!!!\n\n" unless($flag);
			
	$query->finish;
	}
		
	dbDisconnect($dbh) if $close_dbh;
	
	return ($CombsInDataset); # A hash 0f $Hash->{species}{protein_family}=number of times seen.
}


1;
__END__

