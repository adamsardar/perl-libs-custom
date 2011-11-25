package Supfam::Assignments;
require Exporter;

=head1 NAME

Supfam::Assignments.pm

=head1 SYNOPSIS

Functions relating to parsing HMM assignments from SUPERFAMILY (so likely from HMMER)
use Supfam::Assignments;

=head1 AUTHOR

Adam Sardar (Adam.Sardar@bristol.ac.uk)

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 SEE ALSO

=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use DBI;
use Data::Dumper;
use Term::ProgressBar;
use Math::Combinatorics;

use Database::DBfunc;

=pod
=head2 Methods
=over 4
=cut


sub ParseEnsembleAssignments($) {

	my ($AssignmetsFile) = @_;
	
	my $EnsembleAssignmentsHash ={};
	my $EnsembleUnassigned =[];
	
	open FH, "<$AssignmetsFile";
	
	my $dbh = dbConnect();
	my $sth = $dbh->prepare("SELECT sf FROM model WHERE model = ?;");
	
	while (my $line = <FH>){
		
		$line =~ m/^([\w\-]*)\s*([\w\-]*)\s*([\w\-]*)\s*([\w\.\-]*)\s*([\w\-]*)\s*([\w-]*)\s*([\w-]*)\s*([\w-]*)\s*/; 
		#Example lines:
		# ENSP00000386176 0049403 809-842 7.85e-05        9       ESCGLCLKADPDFACGWCqgPGQCTLRQHCPAQE      0.0031  112113  103576
		# ENSP00000386437 -       -       -       -       -       -       -       -
		
		my ($ENSPid,$SFModel,$Region,$Eval,$num,$AssSequence,$Decimal,$num2,$num3) = ($1,$2,$3,$4,$5,$6,$7.$8,$9) #Talk to David and find out what num, num 2 and num3 are, then update variable names
		
		
		unless($SFModel =~ m/\-/){
			
		$EnsembleAssignmentsHash->{$ENSPid} = {} unless(exists($EnsembleAssignmentsHash->{$ENSPid}));
		$EnsembleAssignmentsHash->{$ENSPid}{$Region} = {} unless(exists($EnsembleAssignmentsHash->{$ENSPid}{$Region}));
		
		@{$EnsembleAssignmentsHash->{$ENSPid}{$Region}}{qw(SFModel Eval num1 Sequence Decimal num2 num3)}=($SFModel,$Eval,$num,$AssSequence,$Decimal,$num2,$num3) #Update Hash using a hash slice
		
		$sth->execute($SFModel);
		my $superfamily = $sth->fetchrow_array();
		$EnsembleAssignmentsHash->{$ENSPid}{$Region}{'superfamily'} = $superfamily;
		
		}else{
			
			push(@$EnsembleUnassigned,$ENSPid);
		}
	}
	
	close FH;
	
	dbDisconnect($dbh); 
	
	return ($EnsembleAssignmentsHash,$EnsembleUnassigned);
}

=pod Ensemble2Archs
=head2 Methods
Given a hash created in &ParseEnsembleAssignments, &Ensemble2Archs will produce a hash of ENSPid=>DomArch
=over 4
=cut

sub Ensemble2Archs ($$){
	
	my ($ENSPFullFastaFile,$EnsembleAssignmentsHash)
	
	my $EnsembleFullHash = {}; # This will be a hash of $->{ENSPid}={sequence , Assignments, length, Dom Comb Arch} - maybe consider adding a hash ref of supradomains to this
	
	open ENSPFASTA, "<$ENSPFullFastaFile";
	
	my $FASTAstring; #Because FASTA files are a massive pain in the arse, this string will be the concatenation of multiple lines
	my $ENSPid = 'NULL';
	
	while (my $line = <ENSPFASTA>){
		
		if($line =~ m/^\>(.*)$/){ #Start of a new FASTA entry
			
			my $FullFASTAString = $FASTAstring;
			$EnsembleFullHash->{$ENSPid}={};
			$EnsembleFullHash->{$ENSPid}{'FullFASTA'}=$FullFASTAString;
			$EnsembleFullHash->{$ENSPid}{'FASTAlength'}=length($FullFASTAString);
			
			$ENSPid = $1;			
			$FASTAstring = '';
		
		}else{
			
			chomp($line);
			$FASTAstring .= $line;
		}
	}
	
	delete($EnsembleFullHash->{'NULL'});
	
	$EnsembleFullHash->{$ENSPid}={};
	$EnsembleFullHash->{$ENSPid}{'FullFASTA'}=$FASTAstring;
	$EnsembleFullHash->{$ENSPid}{'FASTAlength'}=length($FASTAstring);

	close ENSPFASTA;

	foreach $ENSPid (keys(%$EnsembleFullHash)){
		
		my $ProteinLength = $EnsembleFullHash->{$ENSPid}{'FASTAlength'};
		my @ProteinRegions = ("1-$ProteinLength"); #This will be populated with assignments
		
		foreach my $AssignedRegion (keys(%{$EnsembleAssignmentsHash->{$ENSPid}})){
			
			$AssignedRegion =~ m/(\d*)\-(\d*)/;
			my ($AssignmenetStart,$AssignmenetEnd) = ($1,$2);
			
			foreach my $FastaPortion (@ProteinRegions){
				
				next if($FastaPortion =~ m/SF:.*/);
				$FastaPortion =~ m/(\d*)\-(\d*)/;
				my ($PortionStart,$PortionEnd) = ($1,$2);
				
				if($AssignmenetStart >= $PortionStart && $AssignmenetStart <= $PortionEnd){
				
					#Delete the region that contains the assingmenet and replace with:
					
					
					#The area before the assignement (unless of length zero)
					my $PreceedingPortion = "$PortionStart".'-'."$AssignmenetStart-1";
					push(@ProteinRegions,$PreceedingPortion) unless($AssignmenetStart == $PortionStart);
					
					#The area after (unless of lenght 0)
					my $PostceedingPortion = "$AssignmenetEnd+1".'-'."$PortionEnd";
					push(@ProteinRegions,$PostceedingPortion) unless($AssignmenetEnd == $PortionEnd);
					
					# The string Modelstaart-SF:(superfamily id)-modelend so as to allow for concatenation of the domain architecture string 
					my $ModelHit = $EnsembleAssignmentsHash->{$ENSPid}{$Region}{'superfamily'};
					my $ModelHitRegionString = "$AssignmenetStart".'-'."SF:$ModelHit".'-'."$AssignmenetEnd";
					
					push(@ProteinRegions,$ModelHitRegionString) unless($AssignmenetEnd == $PortionEnd);
							
					last;	
				}
				#The above block looks fiddly - because it is. Essentially it takes a protein as an unassigned region 1-length(protein). we then cut out of this the portions with assignmenets, leaving an array of
				#gaps (n-m) and assingments (x-SF:ssjdj-y).
				
			}
		}
		
			#sort array.
	
			my @DomainArchitecture; #This will be a concatenated string of domains
			my @gapLengths;
			
			foreach my $ProteinSegment (@ProteinRegions){			
				
				if($ProteinSegment =~ m/^(\d*)\-(\d*)$/){
					
					push(@DomainArchitecture,'_gap_');
					my $gaplength = $2-$1+1;
					push(@gapLengths,$gaplength);
					
				}elsif($ProteinSegment =~ m/^(\d*)\-SF:(\d*)\-(\d*)$/){
					
					push(@DomainArchitecture,$2);
					
				}else{
					die "Script not working as anticipate - non parseable result in ProteinRegions array \n";
				}
				
				my $DomainArchString = join(',',@DomainArchitecture)
				
			}
		
		
		$EnsembleAssignmentsHash->{$ENSPid}{'DomainArchitecture'} = $DomainArchString; #Looking like _gap_,SuperfamilyID,_gap_,SuperfamilyID ...
		$EnsembleAssignmentsHash->{$ENSPid}{'GapLengths'} = \@gapLengths;
		$EnsembleAssignmentsHash->{$ENSPid}{'AssignmentRegions'} = \@ProteinRegions; #This is array looking like ('1-23','24-SF:sfid-134',135-203,'204-SF:sfid-298'...)
	}
	
	
	return($EnsembleAssignmentsHash);
}

#I now have a hash of form: $EnsembleAssignmentsHash->{$ENSPid}..
# {'DomainArchitecture'} = string of domain architecture, {'GapLengths'} = array ref of lengths of the _gap_ assignments,
# {AssignmentRegions} = array ref of which regions in a protein are assigned. The (none _gap_) elements are in turn keys to further entries in the same hash, detialing Evalue etc.

# I assume that I should use the GO_mapping supra table to map GO terms onto proteins. I then need to map SNPs to protein ids and then hand the list to julian.
# Also, think about including Hash's algorithms into all this - is there a deletorious event that might destroy the protein?



1;

