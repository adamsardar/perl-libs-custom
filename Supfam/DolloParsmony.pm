package Supfam::DolloParsmony;
require Exporter;

=head1 NAME

Supfam::DolloParsmony

=head1 SYNOPSIS

A few useful functions for scripting, like easier dumping and reading in of data
use Supfam::DolloParsmony;

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
					PopulateDolloPStateAssingments	
					DolloParsimonyAncestralState
					dolloTraitDecoration
					DOLLOP_Ancestral_Trait_Changes_in_Clade
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use Data::Dumper;
use Supfam::Utils;
use Supfam::TreeFuncsNonBP;

use threads;
use threads::shared;

sub DolloParsimonyAncestralState{ #No prototyping as the number of arguments can be variable
	
	my ($TreeHash,$root,$LeafStatesFile,$nthreads,$TraitLabelsArrayRef) = @_;
	
	my $TraitLabelsProvided = (@_ < 5)?0:1; #Flag for is a list of trait labels provided

	my $TraitLookUp = {};
	
	my @TraitLabels;

	
	open LEAFSTATES, "<$LeafStatesFile";
	
	my $temp = <LEAFSTATES>; #Pass through the first (unimportant) line of phylip traits file detailing the number of traits and number of species
	
	my ($NumberOfSpecies,$NumberOfTraits) = split(/\s+/,$temp);
	
	
	if($TraitLabelsProvided){
	
		@TraitLabels = @$TraitLabelsArrayRef;
		@{$TraitLookUp}{@TraitLabels}=(0 .. $NumberOfTraits-1); #Create a dictionary that maps trait label to position in the string of position - use substr to get the values out
	
	}else{
		
		@TraitLabels = (1 .. $NumberOfTraits);
		@{$TraitLookUp}{@TraitLabels}=(0 .. $NumberOfTraits-1);#Create a dictionary that maps trait label to position in the string of position - use substr to get the values out
	}
	
	die "Different number of species in input file than in the tree! $NumberOfSpecies in file \n\n" unless ($NumberOfSpecies == scalar(@{$TreeHash->{$root}{'Clade_Leaves'}}));
	
	my $NodeNameMapper ={};
	@{$NodeNameMapper}{(map{$TreeHash->{$_}{'node_id'}}(@{$TreeHash->{$root}{'all_Descendents'}}),$root)} = (@{$TreeHash->{$root}{'all_Descendents'}},$root); #Creates a dictionary of nodeID to TreeHash nodeID. These will liekly be the same for leaves, but it's a good idea to allow for flexibility
	
	while (my $line = <LEAFSTATES>){
		
		chomp $line;

		my ($CurrentNode,$AncestralStates) = split(/\s+/,$line);
		
		die "Bad Node Name in genome traits input file: $CurrentNode!\n" unless(exists($NodeNameMapper->{$CurrentNode})); # Bit of error checking

		@TraitLabels = (1 .. length($AncestralStates)) unless ($TraitLabelsProvided); #In case Trait labels was left undefined in input
		
		my $CurrentNodeID = $NodeNameMapper->{$CurrentNode};

		$TreeHash->{$CurrentNodeID}{'DolloPTraitStates'}=$AncestralStates; #A string of binary presecne/abscnece data. Use the 'Trait_String_Poistions_Lookup' entry and substr to extract the data
		$TreeHash->{$CurrentNodeID}{'DolloP_Trait_String_Poistions_Lookup'} = $TraitLookUp; #Have a copy of the trait lookup pointer in every node - it just makes things a little easier and for not much of a memory overhead 
	}
	
close LEAFSTATES;

	##Parsed all the leaf states in. Now to work out all the other information.
	
	my @TreeLeaves = @{$TreeHash->{$root}{'Clade_Leaves'}};
	my @AllNodes = (@{$TreeHash->{$root}{'all_Descendents'}},$root);

	my (undef,undef,undef,$InternalNodes) = IntUnDiff(\@TreeLeaves,\@AllNodes); #Ignore all but final value returned # InUnDiff returns ($Union,$Intersection,$ListAExclusive,$ListBExclusive)

	print STDERR "Parsed Leaf States\n";

	map{$TreeHash->{$_}{'DolloPTraitStates_TempArray'}= [(undef) x $NumberOfTraits]; $TreeHash->{$_}{'DolloP_Trait_String_Poistions_Lookup'} = $TraitLookUp;}@$InternalNodes;
	
	#'DolloPTraitStates_TempArray' is a temporay list to store prescence of a trait. Note that undef is used to save on memory. A trait position will be set as one if the trait exists, left as undef if not
	
	my $TraitDetailsHash = {}; #A hash storing information about which leaves a trait is found in, the MRCA etc.
	
	foreach my $Trait (keys(%{$TreeHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}})){
			
		my $LeavesWithTrait = DolloPLeavesWithTrait($TreeHash,$root,$Trait); #ArrayRef of leaves with given trait.

		my $TraitPosition = $TraitLookUp->{$Trait};
		
		my $MRCA = FindMRCA($TreeHash,$root,$LeavesWithTrait);
		
		$TraitDetailsHash->{$Trait}=[$TraitPosition,$LeavesWithTrait,$MRCA];
	}

	my @Jobs = @AllNodes;
	my $JobHash = {};
	
	my $count=0;
	
	while(my $Job = shift(@Jobs)){
		
		my $Index = $count++ % $nthreads ; 
		$JobHash->{$Index}=[] unless(exists($JobHash->{$Index})); 
		push(@{$JobHash->{$Index}},$Job);
	}

	my @Threads;
	
	foreach my $iterator (keys(%$JobHash)){
		
		my $INPUT = {};
		@{$INPUT}{('TreeHash','TratDetails','JobHash','Iterator')} = ($TreeHash,$TraitDetailsHash,$JobHash,$iterator);

		my $Thr = threads->new(\&Parallel_DolloPStates, $INPUT);
		
		push(@Threads,$Thr);
	}
		
	map{my $ReturnedResults = $_->join(); map{$TreeHash->{$$_[0]}{'DolloPTraitStates'} = $$_[1];}@$ReturnedResults }@Threads;
	
	return(1);
}

=pod
=item *DolloParsimonyAncestralState($TreeHash,$root,$LeafStatesFile,$TraitLabelsArrayRef)

An implementation of the dollo parsimony principle to finding ancestral states. Tries to be effecient with memomory, but no promises (particularly when running many threads).

=cut

sub Parallel_DolloPStates($){
	
	my ($INPUT) = @_;
	
	my ($TreeHash,$TraitDetailsHash,$JobHash,$iterator) = @{$INPUT}{('TreeHash','TratDetails','JobHash','Iterator')} ;
	
	my $Returning = [];
			
	foreach my $Node (@{$JobHash->{$iterator}}){
		
			my $PASS = {}; #PASS I think stands for parallel assignment
			@{$PASS}{('TreeHash','TraitDetails','Node')}=($TreeHash,$TraitDetailsHash,$Node);
		
			my $Reults = PopulateDolloPStateAssingments($PASS);
			push(@$Returning,$Reults);
	}
	
	return $Returning;
}

=pod
=item *Parallel_DolloPStates($TreeHash,$TraitDetailsHash,$JobHash,$iterator)

A quick little wrapper around PopulateDolloPStateAssingments to make parallelisation easier. Please don't call directly.

=cut


sub PopulateDolloPStateAssingments($){
	
	my ($PASS) = @_;

	my ($TreeHash,$TraitDetailsHash,$node) = @{$PASS}{('TreeHash','TraitDetails','Node')};
		
	return([$node,$TreeHash->{$node}{'DolloPTraitStates'}]) if($TreeHash->{$node}{'is_Leaf'}); #Leaves already have their states assigned, so return immediately
	
	my @TempStateArray = (0) x scalar(keys(%$TraitDetailsHash)); #Create an array of zeros that will be updated
	
	foreach my $trait (keys(%$TraitDetailsHash)){
		
		my ($TraitPosition,$LeavesWithTrait,$MRCA) = ($TraitDetailsHash->{$trait}[0],$TraitDetailsHash->{$trait}[1],$TraitDetailsHash->{$trait}[2]);

		my @AllDescendents = @{$TreeHash->{$node}{'all_Descendents'}};
		my @CladeLeaves = @{$TreeHash->{$node}{'Clade_Leaves'}};
		
		next if(grep{m/$MRCA/}@AllDescendents); #If the MRCA of the trait is below this node, then we assume that it must have been created later
		next if (scalar(@$LeavesWithTrait) == 0); #There shoudln't really be any of these if the input file was appropriately formatted, but just in case.
		
		my (undef,$Intersection,undef,undef) = IntUnDiff($LeavesWithTrait,\@CladeLeaves); # InUnDiff returns ($Union,$Intersection,$ListAExclusive,$ListBExclusive)
	
		$TempStateArray[$TraitPosition] = 1 if(scalar(@$Intersection)); #i.e. if any leaves below the node being studied have the trait, then we assume it present here.
	}
	
	my $AncState = join('',@TempStateArray); #Stringlify the states of the node
	
	my $Results = [$node,$AncState];
	
	return $Results; 
}

=pod
=item *PopulateDolloPStateAssingments($TreeHash,$Node,$TraitPosition,$LeavesWithTrait)

Function to assign dollo parsimony states to a single node, $Node. Please don't call directly.

=cut

sub dolloTraitDecoration($$$$$); #Recursive function
sub dolloTraitDecoration($$$$$){
	
	my ($TreeHash,$Node,$FILEHANDLE,$LeavesWithTrait,$FormatString) = @_;
	
	my @NodeLeafDescendents;
	
	unless($TreeHash->{$Node}{'is_Leaf'}){
	
		@NodeLeafDescendents = @{$TreeHash->{$Node}{'Clade_Leaves'}};
	}else{
		
		@NodeLeafDescendents = ($Node);
	}
	
	
	my (undef,$Intersection,undef,undef) = IntUnDiff($LeavesWithTrait,\@NodeLeafDescendents);
		
		
	if(scalar(@$Intersection)){
		
		unless($TreeHash->{$Node}{'is_Leaf'}){
		
			foreach my $Descendent (@{$TreeHash->{$Node}{'each_Descendent'}}){
				
				dolloTraitDecoration($TreeHash,$Descendent,$FILEHANDLE,$LeavesWithTrait,$FormatString);
			}
		}else{
			
			return(1);
			
		}
			
	}else{
			
		my $TraitAbsentCladeString = join(' ',map{$TreeHash->{$_}{'node_id'}}@NodeLeafDescendents);
	
		print $FILEHANDLE $FormatString.$TraitAbsentCladeString."\n";
			
		return(1);
	}
	
	
}

=pod
=item *dolloTraitDecoration($TreeHash,$MRCA,FILEHANDLE,$LeavesWithTrait,$FormatString)

Function to assign dollo parsimony estimate of a single trait on a given tree and output
the clade's recognised to a file FILEHANDLe using the style given in $FormatString.

=cut


=item * DOLLOP_Ancestral_Trait_Changes_in_Clade($TreeCacheHash,$node,$traitlabelsarray)

After running DolloParsimonyAncestralState to calculate dollo parsimony ancestral states, this function will calculate the number of traits that have been created or deleted along each branch.
These are stored in:

	$TreeCacheHash->{$node}{'DOLLOP_Number_Created'}
	$TreeCacheHash->{$node}{'DOLLOP_Number_Deleated'}

Also present for each node are:

	$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Created'}
	$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Deleated'}

=cut

sub DOLLOP_Ancestral_Trait_Changes_in_Clade{ #Left out the prototyping arguments as its a recursive function. Be careful
	
	my ($TreeCacheHash,$node,$root,$traitlabelsarray) = @_;
	
	my $TraitLabelsProvided = (@_ < 4)?0:1; #Flag for is a list of trait labels provided - these are the traits that we wish to act upon. If nothign is provided, then all traits will be used
	@$traitlabelsarray = (keys(%{$TreeCacheHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}})) unless($TraitLabelsProvided); #If no traits were provided, assume that all labels are wanted
	
	die "No entries in '$TreeCacheHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}'. Run DolloParsimonyAncestralState first before using this function." unless(exists($TreeCacheHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}));
	
	my @AncestorStates;
	my $TraitsPositionsLookup = $TreeCacheHash->{$node}{'DolloP_Trait_String_Poistions_Lookup'}; #All traits are stored as a logn string like 01001010010. We use substr to extract the one we want. This dictionary simply maps from trait to position
	
	if($node eq $root){
		
		@AncestorStates = (0) x scalar(@$traitlabelsarray); #If the node is the root, then all prescent traits are assumed to have been created beforehand. So it's ancestor will be all 0's.
	
	}else{
	
		my $Ancestor = $TreeCacheHash->{$node}{'ancestor'};
		my $AncestorStateString = $TreeCacheHash->{$Ancestor}{'DolloPTraitStates'};
		@AncestorStates = map{substr($AncestorStateString,$_,1)}(@{$TraitsPositionsLookup}{@$traitlabelsarray});
	}
	
	my $NodeStateString = $TreeCacheHash->{$node}{'DolloPTraitStates'};
	my @NodeStates = map{substr($NodeStateString,$_,1)}(@{$TraitsPositionsLookup}{@$traitlabelsarray});
		
	my ($no_creations,$no_deletions) = (0,0); 
	
	foreach my $TraitIndex (0 .. (scalar(@$traitlabelsarray)-1)){
		
		my ($AncestorState,$NodeState) = ($AncestorStates[$TraitIndex],$NodeStates[$TraitIndex]);
		
		next if ($AncestorState ~~ '?' || $NodeState ~~ '?');
		
		die "Only binary states supported  $AncestorState $NodeState \n" unless ( scalar(grep{$_ == $AncestorState}(0,1)) && scalar(grep{$_ == $NodeState}(0,1)));
	
		if($AncestorState - $NodeState < 0){
	
			$no_creations ++;
			
		}elsif($AncestorState - $NodeState > 0){
	
			$no_deletions ++;
		}
	}
				
	$TreeCacheHash->{$node}{'DOLLOP_Number_Created'} = $no_creations;
	$TreeCacheHash->{$node}{'DOLLOP_Number_Deleted'} = $no_deletions;
	#This is the number of traits created/deleted along the ancestral branch ABOVE the current node 
	
	unless ($TreeCacheHash->{$node}{'is_Leaf'}){
			
		my @Children = @{$TreeCacheHash->{$node}{'each_Descendent'}};
			
		$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Created'} = $no_creations;
		$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Deleted'} = $no_deletions;	
			
		foreach my $Child (@Children){
			
			 my ($Total_no_creations,$Total_no_deletions) = DOLLOP_Ancestral_Trait_Changes_in_Clade($TreeCacheHash,$Child,$root,$traitlabelsarray);
			#Notice the replacement of $node with $Child.
			#Total_no_X is the number of deletions/ creations in branches BELOW this node. 
			#This will be added to the number of deletions along this nodes ancestral branch to give the total in the clade BELOW the branch beneath the ancestor of the node under study

			$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Created'} += $Total_no_creations;
			$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Deleted'} += $Total_no_deletions;
		}
		
		my ($Total_no_creations,$Total_no_deletions) = ($TreeCacheHash->{$node}{'DOLLOP_Total_Number_Created'},$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Deleted'});
		
		return($Total_no_creations,$Total_no_deletions);
			
	}else{
			
		$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Created'} = $no_creations;
		$TreeCacheHash->{$node}{'DOLLOP_Total_Number_Deleted'} = $no_deletions;	
			
		return($no_creations,$no_deletions);
	}
}




1;

