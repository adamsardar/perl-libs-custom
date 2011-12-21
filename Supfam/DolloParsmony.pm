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
	
	my $temp = <LEAFSTATES>; #Pass through the first (unimportant) line of phylip file detailing the number of traits and number of species
	
	my ($NumberOfSpecies,$NumberOfTraits) = split(/\s+/,$temp);
	
	#my $RAxMLTraitsPresent = 1;
	
	#map{$RAxMLTraitsPresent = 0 unless(exists($TreeHash->{$_}{'RAxML_AncestralStates'}));}@{$TreeHash->{$root}{'Clade_Leaves'}}; #No point reparsing a file if the values already exist in memory
	
#unless($RAxMLTraitsPresent){
	
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
		
		die "Bad Node Name in input file: $CurrentNode!\n" unless(exists($NodeNameMapper->{$CurrentNode})); # Bit of error checking

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

	print "Parsed Leaf States";

	map{$TreeHash->{$_}{'DolloPTraitStates_TempArray'}= [(undef) x $NumberOfTraits]; $TreeHash->{$_}{'DolloP_Trait_String_Poistions_Lookup'} = $TraitLookUp;}@$InternalNodes;
	
	#'DolloPTraitStates_TempArray' is a temporay list to store prescence of a trait. Note that undef is used to save on memory. A trait position will be set as one if the trait exists, left as undef if not
	
	my $TraitDetailsHash = {}; #A hash storing information about which leaves a trait is found in, the MRCA etc.
	
	foreach my $Trait (keys(%{$TreeHash->{$root}{'DolloP_Trait_String_Poistions_Lookup'}})){
			
		my $LeavesWithTrait = LeavesWithTrait($TreeHash,$root,$Trait); #ArrayRef of leaves with given trait.

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
	
	print "\n\nCalculated Dollo Parsimony estimates of each trait\n\n";
	
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
		
		next if(grep{m/$MRCA/}@AllDescendents); #If the MRCA of the trait is below this node, then we assume that it must have been created later
		next if (scalar(@$LeavesWithTrait) == 0); #There shoudln't really be any of these if the input file was appropriately formatted, but just in case.
		
		my (undef,$Intersection,undef,undef) = IntUnDiff($LeavesWithTrait,\@AllDescendents); # InUnDiff returns ($Union,$Intersection,$ListAExclusive,$ListBExclusive)
	
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





1;

