package Supfam::Phylip_Ancestral_States_Parser;
require Exporter;

=head1 NAME

Supfam::Phylip_Ancestral_States_Parser

=head1 SYNOPSIS

A very small modules for the parsing of Phylip format ancestral state reconstruction files, whether DOLLOP or MIX
use Supfam::Phylip_Ancestral_States_Parser;

=head1 AUTHOR

Adam Sardar (adam.sardar@bris.ac.uk)

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 SEE ALSO

Supfam::Config.pm

=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
					ParsePhylipAnces
					MapPhylipTree2BioPerl
					MapPhylipTree2TreeCache
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use Supfam::Utils;

sub ParsePhylipAnces($) {
	
	my ($InputFile) = @_;
	
	open my $FH,  "< $InputFile" or die "$? $!";
	
	my $TreeBranchHash = {};
	
	my $AncesStatesFlag = 0; #Flag to check if we're in the correct portion of the file.
	
	my ($BranchChild,$BranchParent,$BranchArrayRef) = ('NULL',undef,undef);
	
	foreach my $line (<$FH>){
		
		next if ($line =~ m/^$/ || $line =~ m/^#/); #Weed out comments and empty lines
		chomp($line); #Remove trailing new line
		
		if($AncesStatesFlag){
			
			my @SplitLine; # This will be the input line split on whitespace
						
			if($line =~ m/\s+(yes|no|maybe|root|\d{1,3})\s+/i){ # Used to catch the comments in the column corresponding to 'Any Steps?' - shows that we are dealing with a new data entry
				
				@SplitLine = split(/\s+/,$line); # example line 'root      1         yes    .1..1 ..1.1'
								
				{no warnings; while($SplitLine[0] =~ m/^$/){ shift(@SplitLine);}} #That phylip file is such a pain in the arse to parse! This rips off the leading whitspaces
				#The warnigns are turned off above so as not to kick out errors from the while condition
				
				$TreeBranchHash->{$BranchChild}[1] = join(@{$TreeBranchHash->{$BranchChild}[1]}) if(exists($TreeBranchHash->{$BranchChild})); #Stringify the array of traits to lower peak memory footprint
				
				$BranchParent = shift(@SplitLine);
				$BranchChild = shift(@SplitLine);
				shift(@SplitLine); #Eliminate the yes/no entry
				
				$BranchChild =~ s/Taxon(\w{2,3})/$1/ if($BranchChild =~ m/Taxon\w{2,3}/);
				
				#Data structure up to this point is: $TreeBranchHash->{$BranchChild}[BranchAncestor,[list of trait changes]] . In order to keep memory usage down, this array is stringified and replaces the array ref
				#So instead the data struct becomes $TreeBranchHash->{$BranchChild}[BranchAncestor,string of trait changes]
				
				$TreeBranchHash->{$BranchChild} = [$BranchParent,[]]; #Data structure is of the form $TreeBranchHash->{$BranchChild}[BranchAncestor,[list of trait changes]] 
				
			}else{
			
				@SplitLine = split(/\s+/,$line); # example line '.1..1 ..1.1'

				{no warnings; while ($SplitLine[0] =~ m/^$/){ shift(@SplitLine);}} #That phylip file is such a pain in the arse to parse! This rips off the leading whitspaces
				#The warnigns are turned off above so as not to kick out errors from the while condition
			}
			
			while (my $StateChanges = shift(@SplitLine)){
									
				$BranchArrayRef = $TreeBranchHash->{$BranchChild}[1]; # $TreeBranchHash->{$BranchChild}[1] refers to the branch connecting $BranchChild to its parent
				#The array that items are pushed onto are the changes along that branch (1 -> creataion, 0 -> deletion, .->no change)
				
				push(@$BranchArrayRef,split('',$StateChanges)); # The array @$BranchChild will look like (0,0,0,.,.,.,,1,1, etc.)
			}
			
		}
		
		$AncesStatesFlag = 1 if($line =~ m/means same as in the node below it on tree/g); #I've chosen an arbitaty piece of comment that wont be repeated elsewhere to initiate the file parsing point
	}
	
	$TreeBranchHash->{$BranchChild}[1] = join(@{$TreeBranchHash->{$BranchChild}[1]}); #So much of this script is horrible because parsing the data file is just plain painful!
	
	
	
	# Returned is:$TreeBranchHash->{$BranchChild}[BranchAncestor,string of trait changes]  (1 -> creataion, 0 -> deletion, .->no change)
	
	close $FH;
	return($TreeBranchHash);
}

sub MapPhylipTree2BioPerl($$$){
	
	my ($AncestralChangesHash,$TreeCacheHash,$root) = @_ ;
		
	my $Bioperl2PhylipMap = {};
	my $Phylip2BioperlMap = {};
		
	my @BPCladeLeaves = @{$TreeCacheHash->{$root}{'Clade_Leaves'}}; #List of leaves as bioperl nodeIDs
		
	$Bioperl2PhylipMap->{$root}='1';
	$Phylip2BioperlMap->{'1'}=$root;
	#Initialise the easy value in the tree mapping root to root
	
	foreach my $BPLeaf (@BPCladeLeaves){
		
		my $LeafName = $TreeCacheHash->{$BPLeaf}{'node_id'};
		my $PhylipLeaf = "$LeafName";
		print $PhylipLeaf."\n";
		
		$Bioperl2PhylipMap->{$BPLeaf}=$PhylipLeaf;
		$Phylip2BioperlMap->{$PhylipLeaf}=$BPLeaf;
		
		my $BioperlAncestor = $TreeCacheHash->{$BPLeaf}{'ancestor'};
		my $PhylipAncestor = $AncestralChangesHash->{$PhylipLeaf}[0];
		
		print "out $BioperlAncestor <- $BPLeaf  BP  $PhylipAncestor <-$PhylipLeaf Phy \n";
		while (! exists($Bioperl2PhylipMap->{$BioperlAncestor})){ #This will at most go up to the root of the tree sicne we set this earlier
			
					print "$PhylipAncestor in\n";

					$Bioperl2PhylipMap->{$BioperlAncestor}=$PhylipAncestor;
					$Phylip2BioperlMap->{$PhylipAncestor}=$BioperlAncestor;
					
					my ($BPnode,$PhyNode) = ($BioperlAncestor,$PhylipAncestor);
					
					$BioperlAncestor = $TreeCacheHash->{$BioperlAncestor}{'ancestor'};
					$PhylipAncestor = $AncestralChangesHash->{$PhylipAncestor}[0];
					print "in $BioperlAncestor <- $BPnode  BP  $PhylipAncestor <-$PhyNode Phy \n";
		}
	}
	
	return($Bioperl2PhylipMap,$Phylip2BioperlMap);
}

sub MapPhylipTree2TreeCache{
	
	my ($TreeHash,$PhylipAncestralHash,$TreeHash2PhylipDictionary,$Phylip2TreeHashDictionary,$TreeHashnode) = @_;
	
	my $PhylipAncestor;
   my $PhylipNode;
	
	unless($TreeHash->{$TreeHashnode}{'is_Leaf'}){
		
		my @PhylipNodeLabel; #This will be what the descendent nodes think that his node is in the phylip schema
		
		foreach my $TreeHashdescendent (@{$TreeHash->{$TreeHashnode}{'each_Descendent'}}){
	
			$PhylipNode = MapPhylipTree2TreeCache($TreeHash,$PhylipAncestralHash,$TreeHash2PhylipDictionary,$Phylip2TreeHashDictionary,$TreeHashdescendent);
			
			push(@PhylipNodeLabel,$PhylipNode);
			
			print $TreeHashdescendent.scalar(@{$TreeHash->{$TreeHashnode}{'each_Descendent'}})."\n";
		
			
		}
		
		my $Number = grep{m/^$PhylipNodeLabel[0]$/}@PhylipNodeLabel;

		print join("\n",@PhylipNodeLabel);
		print "\n";
		
		die "Disagreement as to what this node is in phylip - error somewhere! $Number  $TreeHashnode <- TreeHash  $PhylipNode <- Phylip " unless ((grep{m/^$PhylipNodeLabel[0]$/}@PhylipNodeLabel) == scalar(@PhylipNodeLabel));
		
		$TreeHash2PhylipDictionary->{$TreeHashnode} = $PhylipNode;
		$Phylip2TreeHashDictionary->{$PhylipNode} = $TreeHashnode;
		
		$PhylipAncestor = $PhylipAncestralHash->{$PhylipNode}[0];
		
	}else{
				
		print 'Leaf'."\n";
		#Node is a leaf, for which we know that phylip will have a direct 1-1 mapping
		
		die "Leaves don't appear to agree!" unless(exists($PhylipAncestralHash->{$TreeHashnode}));
				
		$PhylipNode = $TreeHashnode;
				
		$TreeHash2PhylipDictionary->{$TreeHashnode} = $PhylipNode;
		$Phylip2TreeHashDictionary->{$PhylipNode} = $TreeHashnode;
		
		$PhylipAncestor = $PhylipAncestralHash->{$PhylipNode}[0];
	}
	
	return($PhylipAncestor); #For the very first call, this should be root.
	
}




1;