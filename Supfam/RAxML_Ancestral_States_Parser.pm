package Supfam::RAxML_Ancestral_States_Parser;
require Exporter;

=head1 NAME

Supfam::RAxML_Ancestral_States_Parser

=head1 SYNOPSIS

A few useful functions for scripting, like easier dumping and reading in of data
use Supfam::RAxML_Ancestral_States_Parser;

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
					RAxMLAncestralMarginalProbabilities_FileParser
					RAxMLAncestralMarginalStates_FileParser
					DeletionsTreeByMLProbabilities
					RAxML_Ancestral_Trait_Changes_in_Clade
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use Supfam::TreeFuncs;
use Supfam::Utils;

use Bio::Tree::TreeI;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;

sub RAxMLAncestralMarginalProbabilities_FileParser{
	
	my ($RAXMLAncestralMarginalProbabilitiesFile,$SpeciesTraitsFile,$TreeCacheHash,$TraitLabelsArrayRef) = @_;
	
	my $TraitLabelsProvided = (@_ < 4)?0:1; #Flag for is a list of trait labels provided
	
	my $TraitLookUp = {}; #This will be a way of finding out where a trait in a strong is stored
	my @TraitLabels;
	
	my @NodeNames = map{$TreeCacheHash->{$_}{'node_id'}}(keys(%$TreeCacheHash));
	
	my $NodeNameHash = {};
	@{$NodeNameHash}{@NodeNames}=undef;
	#Creates a hash for fast lookups - useful when we loop through the file and identify lines that correspond to particular entries

	die "Non-unique node IDs. Bad plan!" unless (scalar(keys%$NodeNameHash) == scalar(@NodeNames));
	
	open SPECIESSTATES, "<$SpeciesTraitsFile" or die $!.$?;
	
	my $temp = <SPECIESSTATES>; #skip past a useless file header
	
	my ($NumberOfSpecies,$NumberOfTraits) = split(/\s+/,$temp);
	my $NoOfStates = 2;
	
	print "Input Number of traits = $NumberOfTraits\n";
	print "Input Number of species = $NumberOfSpecies\n";
	print "Input Number of states = $NoOfStates\n";
	
	if($TraitLabelsProvided){
		
		die "There should be the same number of trait labels as traits!\n" unless (scalar(@$TraitLabelsArrayRef) == $NumberOfTraits);# Quick bit of error checking
		
		@TraitLabels = @$TraitLabelsArrayRef;
		map{$TraitLookUp->{$TraitLabels[$_]} = [$_*8,$_*8+7] }(0 .. $NumberOfTraits-1); #For each trait, have an array of trait start and trait stop in string.
		#$TraitLookUp is a dictionary that maps trait label to position in the string of position - use substr to get the values out
	}else{
		
		@TraitLabels = (1 .. $NumberOfTraits);
		map{$TraitLookUp->{$TraitLabels[$_]} = [$_*8,$_*8+7] }(0 .. $NumberOfTraits-1); #For each trait, have an array of trait start and trait stop in string.
		#$TraitLookUp is a dictionary that maps trait label to position in the string of position - use substr to get the values out
		#$TraitLookUp->{trait_label} = [start end]
	}

	#BioPerlId->Node Name to internal node id (assigned when $TreeCacheHash was created)
	my $NodeName2BPNodeIDMapper = {};
	map{$NodeName2BPNodeIDMapper->{$TreeCacheHash->{$_}{'node_id'}}=$_}(keys(%$TreeCacheHash));
	
	my $CurrentNode; #Throughout the loop this will be the node under analysis. Suggest a better format to Alexis of RAxML
	my $AncestralState; #This will be an array ref concatenated string of multistate probabilities. So length 2 in the binary case
 	
 	my ($State0,$State1); #ASSUMPTION! At current RAxML only deals with binary state ancestral reconstruction. If this changes, then this part of the script will need altering
 		
	while (my $line = <SPECIESSTATES>){
		
		next if($line =~ m/^#/ || $line =~ m/^$/); #Trim out blank lines and comments
		
		chomp($line);
		
		my ($species,$traits)= split(/\s+/,$line);
		my $BPCurrentNode = $NodeName2BPNodeIDMapper->{$species};
		
		($State0,$State1) = ([(("0.000000")x$NoOfStates)],[(("0.000000")x$NoOfStates)]);
		
		my $TraitCount = 0;	
		map{if($_){$State0->[$TraitCount] = "1.000000";}else{$State0->[$TraitCount] = "1.000000";} $TraitCount++; }(split('',$traits));
			
		my $traitsarray = [(join('',@$State0),join('',@$State0))];

		$TreeCacheHash->{$BPCurrentNode}{'RAxML_AncestralProbabilities'} = $traitsarray;
		$TreeCacheHash->{$BPCurrentNode}{'RAxML_Trait_String_Poistions_Lookup'} = $TraitLookUp; #Have a copy of the trait lookup pointer in every node - it just makes things a little easier and for not much of a memory overhead 
	}
	
	close SPECIESSTATES;
	

	open RAXMLSTATES, "<$RAXMLAncestralMarginalProbabilitiesFile" or die $!.$?;
	
	while (my $line = <RAXMLSTATES>){
		
		next if($line =~ m/^#/ || $line =~ m/^\s*$/); #Trim out blank lines and comments
		
		no warnings; #Turns off warnings for the sake of the exists check in the following if() statement
			
		chomp($line);

			if (exists($NodeNameHash->{$line})){
				
				$CurrentNode=$line; #Switch the node under anlaysis.			
				
				if($AncestralState !~ undef){
					
					my $State0String = join('',@$State0);
					my $State1String = join('',@$State1);
					
					@$AncestralState = ($State0String,$State1String)
				}
				
				($State0,$State1) = ([],[]);
				$AncestralState = [];
				my $BPCurrentNode = $NodeName2BPNodeIDMapper->{$CurrentNode};
				$TreeCacheHash->{$BPCurrentNode}{'RAxML_AncestralProbabilities'} = $AncestralState;
				$TreeCacheHash->{$BPCurrentNode}{'RAxML_Trait_String_Poistions_Lookup'} = $TraitLookUp; #Have a copy of the trait lookup pointer in every node - it just makes things a little easier and for not much of a memory overhead 
				
			}else{
				
				my @MultiStateProbablities = split(/\s+/,$line);
				die "Script only handles binary state reconstruction at line $." if(scalar(@MultiStateProbablities) != 2);	
				#Change the states to shorter ints (6d.p.)
								
				map{unless($_ =~ m/(0|1)\.\d{6}/){ die "Expecting marginal probabilities to be of the for 0.123456 with six digits agter dp! Got $_\n"}}@MultiStateProbablities;
				#Expected 8 carachters per probability
				
				#push onto states
				push(@$State0,$MultiStateProbablities[0]);
				push(@$State1,$MultiStateProbablities[1]); 
			}
	}
	
	
	my $State0String = join('',@$State0);
	my $State1String = join('',@$State1);
	@$AncestralState = ($State0String,$State1String);
	#Process the last trait in the file
	
	my $RootBP = $NodeName2BPNodeIDMapper->{'ROOT'};
	
	my $MeasuredNoOfTraits = length($TreeCacheHash->{$RootBP}{'RAxML_AncestralProbabilities'}[0])/8;

	open FH, "> file.traits" or die $!;
	
	print FH $TreeCacheHash->{$RootBP}{'RAxML_AncestralProbabilities'}[0];
	
	close FH;
	
	#Using the root to estimate the number of states and the number of traits should be fine, but we'll use the first level down just to be on the safe side

	print "Number of states = $NoOfStates\n";
	print "Number of traits = $MeasuredNoOfTraits\n";
	
	
	my @RootDescendents = @{$TreeCacheHash->{$RootBP}{'each_Descendent'}};
	
	close RAXMLSTATES;

	print "Finally!\n";
	
	return(1);
	
	#A script to take a list of genomes (strains or whatever) and output some summary data regarding the domain architectures in common.
	#ally updates $TreeChacehHash. YEs this will use alot of memory. But, if run over the whole of SUPERFAMILY on
	#complete tree of life this will take up 2GB. Not so bad. The other option is to migrate to an object based system which only pulls in the
	#relevant stuff on the fly. Or a database strucutre ...
}

=pod
=item * RAxMLAncestralMarginalStates_FileParser(RAXMLAncestralMarginalProbabilitiesFileHandle,$TreeCacheHash)

Using the predifned data structure $TreeChacheHash from 'pacakge Supfam::TreeFuncs', we create a new entry:


$TreeChacheHash->{BioPerlNodeID}{'RAxML_AncestralProbabilities'}=[[p(state0 for trait 1),p(state1 for trait 1),..],[p(state0 for trait 2),p(state1 for trait 2),..],...,[p(state0 for trait n),p(state1 for trait n),..]]

This uses alot of memory. But for all of SUPERFAMILY, under binary states, for all include='y' genomes, this is 2GB of disk. Which is just about OK.
Will consider an OO method that create ancestral states on the fly if this is too memory intensive.

=cut

sub RAxMLAncestralMarginalStates_FileParser{

	my ($RAXMLAncestralMarginalStatesFile,$TreeCacheHash,$LeafTraitsFile,$TraitLabelsArrayRef) = @_;
		
	my $TraitLabelsProvided = (@_ < 4)?0:1; #Flag for is a list of trait labels provided
	
	my $TraitLookUp = {};
	my @TraitLabels;
	@TraitLabels = @$TraitLabelsArrayRef if ($TraitLabelsProvided);
	my $NumberOfTraits = undef;
	@{$TraitLookUp}{@TraitLabels}=(0 .. scalar(@TraitLabels)-1) if ($TraitLabelsProvided); #Create a dictionary that maps trait label to position in the string of position - use substr to get the values out
	$NumberOfTraits = scalar(keys(%$TraitLookUp)) if($TraitLabelsProvided); # I'm placing a fair amoung of trust in RAxML that the number of traits is constant throughout the file, but this is a possible place for a bug
	
	my @NodeNames = map{$TreeCacheHash->{$_}{'node_id'}}(keys(%$TreeCacheHash));
	my $NodeNameHash = {};@{$NodeNameHash}{@NodeNames}=undef;
	#Creates a hash for fast lookups - useful when we loop through the file and identify lines that correspond to known nodes it the tree (stops us having to loop through a list instead)
	die "Non-unique node IDs. Bad plan!" unless (scalar(keys%$NodeNameHash) == scalar(@NodeNames)); #Error check
	
	#BioPerlId->Node Name to internal node id (assigned when $TreeCacheHash was created)
	my $BPNodeID2NodeNameMapper = {}; map{$BPNodeID2NodeNameMapper->{$TreeCacheHash->{$_}{'node_id'}}=$_}(keys(%$TreeCacheHash));
	
	open RAXMLSTATES, "<$RAXMLAncestralMarginalStatesFile";
	
	while (my $line = <RAXMLSTATES>){
		
		chomp $line;

		my ($CurrentNode,$AncestralStates) = split(/\s+/,$line);
		
		if ($. == 1 && ! $TraitLabelsProvided){
			
			@{$TraitLookUp}{(1 .. length($AncestralStates))}=(0 .. length($AncestralStates)-1); #Create a dictionary that maps trait label to position in the string of position - use substr to get the values out		
			#As this only gets created if no traitlabels were specified, then it's a simple index from 1 through lenght(state line in phylip)
			$NumberOfTraits = scalar(keys(%$TraitLookUp));
		}
		
		
		die "Bad NodeID in input file!\n" unless(exists($NodeNameHash->{$CurrentNode})); # Bit of error checking
		
		my $BPCurrentNodeID = $BPNodeID2NodeNameMapper->{$CurrentNode};

		$TreeCacheHash->{$BPCurrentNodeID}{'RAxML_AncestralStates'}= $AncestralStates;
		$TreeCacheHash->{$BPCurrentNodeID}{'Trait_String_Poistions_Lookup'} = $TraitLookUp; #Have a copy of the trait lookup pointer in every node - it just makes things a little easier and for not much of a memory overhead 
	}
		
	close RAXMLSTATES;
	
	print STDERR "Parsed internal states\n";
	
	open LEAFSTATES, "<$LeafTraitsFile";
	
	my $temp = <LEAFSTATES>; #Pass through the first (unimportant) line of phylip file
	
	while (my $line = <LEAFSTATES>){
		
		chomp $line;

		my ($CurrentNode,$AncestralStates) = split(/\s+/,$line);
		
		die "Bad NodeID in input file!\n" unless(exists($NodeNameHash->{$CurrentNode})); # Bit of error checking

		@TraitLabels = (1 .. length($AncestralStates)) unless ($TraitLabelsProvided); #In case Trait labels was left undefined in input
		
		my $BPCurrentNodeID = $BPNodeID2NodeNameMapper->{$CurrentNode};

		$TreeCacheHash->{$BPCurrentNodeID}{'RAxML_AncestralStates'}=$AncestralStates; #A string of binary presecne/abscnece data. Use the 'Trait_String_Poistions_Lookup' entry and substr to extract the data
		$TreeCacheHash->{$BPCurrentNodeID}{'Trait_String_Poistions_Lookup'} = $TraitLookUp; #Have a copy of the trait lookup pointer in every node - it just makes things a little easier and for not much of a memory overhead 
	}
		
	close LEAFSTATES;
	
	print STDERR "Parsed leaf states\n";
	
	return($TraitLookUp); #Function updates $TreeCacheHash so that one can read off the ancestral states
}


=pod
=item * RAxMLAncestralMarginalStates_FileParser()

Using the predifned data structure $TreeChacheHash from 'pacakge Supfam::TreeFuncs', we create a new entry:


$TreeChacheHash->{BioPerlNodeID}{'RAxML_AncestralStates'}={'trait0'='state','trait1'=>state1}

If a list of trait names is provided (of length equal to the number of states in the RAxML output), then these will be used as labels. Else 1 through n shall be used (this allows for easy comparison of traits)

=cut

sub DeletionsTreeByMLProbabilities{
	
	my ($TreeCacheHash,$root)= @_;
	
	#$tree is the bioperl object of the newick tree. $TreechacheHash is a hash full of additional information
	
	my @BPNodeIDs = keys(%$TreeCacheHash);
	my @LabelledNodeIDs = map{$TreeCacheHash->{$_}{'node_id'}}@BPNodeIDs;

	map{die "Run &RAxMLAncestralMarginalProbabilities_FileParser() first on TreeCacheHash in order to create ancestral state entries\n" unless(exists($TreeCacheHash->{$_}{'RAxML_AncestralProbabilities'}))}@BPNodeIDs;
	#Quick bit of error checking
	
	my @NonRootNondes = grep{!/$root/}@BPNodeIDs;
	
	my $RandomNode = $NonRootNondes[int(rand(@NonRootNondes))];
	
	my $NumberOfStates = scalar(@{$TreeCacheHash->{$RandomNode}{'RAxML_AncestralProbabilities'}[0]});
	my $NumberOfTraits = scalar(@{$TreeCacheHash->{$RandomNode}{'RAxML_AncestralProbabilities'}});

	print "Number of states = $NumberOfStates\n";
	print "Number of traits = $NumberOfTraits\n";

	my $ArrayPoint = $TreeCacheHash->{$RandomNode}{'RAxML_AncestralProbabilities'};
	EasyDump("Dump", $ArrayPoint);
	
	die "At current, only a binary state representation is supported" unless($NumberOfStates == 2);
	#$TreeCacheHash->{node}{'AncestralProbabilities'} is a 2D array:[[p(state0 for trait 1),p(state1 for trait 1),..],[p(state0 for trait 2),p(state1 for trait 2),..],...,[p(state0 for trait n),p(state1 for trait n),..]]
	
	foreach my $node (keys(%$TreeCacheHash)){
		
		next if($node eq $root);
		#Root has no ancestor by definition, so we should leave it out of analysis
		
		my $NodeParent = $TreeCacheHash->{$node}{'ancestor'};

		my @nodestate;
		my @ancestorstate;
		map{push(@nodestate,$TreeCacheHash->{$node}{'RAxML_AncestralProbabilities'}[$_][1])}(0 .. $NumberOfTraits-1);
		map{push(@ancestorstate,$TreeCacheHash->{$NodeParent}{'RAxML_AncestralProbabilities'}[$_][1])}(0 .. $NumberOfTraits-1);
		#We are only studying the present state - so we look at [trait][1]		
		
		my @StateChange = map{$ancestorstate[$_] - $nodestate[$_]}(0 .. $NumberOfTraits-1);
		#This is the change in marginal probabilities, positive and negative, over the branch, per 
		
		my @DeletionsOverBranch = grep{$_ >= 0}@StateChange;
		#We are only interested in deletions along the branch - i.e. a higher probability of seeing the trait in an ancestor than the descendent 
		
		my $NormalisedBranch = 0; map{$NormalisedBranch += $_}@DeletionsOverBranch;
		$NormalisedBranch=$NormalisedBranch/$NumberOfTraits; #normalise so that we don't have a huge number as a branch length
	
		$TreeCacheHash->{$node}{'branch_length'} = $NormalisedBranch;
	}
	
	my $NewickTree = TreeHash2Newick($TreeCacheHash,$root);

	return($NewickTree);
}

=pod
=item * DeletionsTreeByMLProbabilities($TreeCacheHash,$root,$tree)

Using the predifned data structure $TreeChacheHash from 'pacakge Supfam::TreeFuncs', 

we calculate the probability of deletion along and edge and sum this for each trait. This is then proportional to the branch length in returned tree returned (same topology as before)

If a list of trait names is provided (of length equal to the number of states in the RAxML output), then these will be used as labels. Else 1 through n shall be used (this allows for easy comparison of traits)

=cut



=item * RAxML_Ancestral_Trait_Changes_in_Clade($TreeCacheHash,$node,$traitlabelsarray)

After running RAxMLAncestralMarginalStates_FileParser on a RAxML ancestral output, this function will calculate the number of rounded traits (i.e thresholded probabilities) that have been created or deleted along each branch.
These are stored in:

	$TreeCacheHash->{$node}{'RAxML_Number_Created'}
	$TreeCacheHash->{$node}{'RAxML_Number_Deleated'}

Also present for each node are:

	$TreeCacheHash->{$node}{'RAxML_Total_Number_Created'}
	$TreeCacheHash->{$node}{'RAxML_Total_Number_Deleated'}	

=cut


sub RAxML_Ancestral_Trait_Changes_in_Clade{ #Left out the prototyping arguments as its a recursive function. Be careful
	
	my ($TreeCacheHash,$node,$root,$traitlabelsarray) = @_;
	
	my $TraitLabelsProvided = (@_ < 4)?0:1; #Flag for is a list of trait labels provided - these are the traits that we wish to act upon. If nothign is provided, then all traits will be used
	@$traitlabelsarray = (keys(%{$TreeCacheHash->{$root}{'Trait_String_Poistions_Lookup'}})) unless($TraitLabelsProvided); #If no traits were provided, assume that all labels are wanted
		
	die "No entries in '$TreeCacheHash->{$root}{'Trait_String_Poistions_Lookup'}'. Run RAxMLAncestralMarginalStates_FileParser first before using this function." unless(exists($TreeCacheHash->{$root}{'Trait_String_Poistions_Lookup'}));
	
	my @AncestorStates;
	my $TraitsPositionsLookup = $TreeCacheHash->{$node}{'Trait_String_Poistions_Lookup'}; #All traits are stored as a logn string like 01001010010. We use substr to extract the one we want. This dictionary simply maps from trait to position
	
	if($node eq $root){
		
		@AncestorStates = (0) x scalar(@$traitlabelsarray); #If the node is the root, then all prescent traits are assumed to have been created beforehand. So it's ancestor will be all 0's.
	}else{
	
		my $Ancestor = $TreeCacheHash->{$node}{'ancestor'};
		my $AncestorStateString = $TreeCacheHash->{$Ancestor}{'RAxML_AncestralStates'};
		@AncestorStates = map{substr($AncestorStateString,$_,1)}(@{$TraitsPositionsLookup}{@$traitlabelsarray});
	}
	
	my $NodeStateString = $TreeCacheHash->{$node}{'RAxML_AncestralStates'};
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
				
	$TreeCacheHash->{$node}{'RAxML_Number_Created'} = $no_creations;
	$TreeCacheHash->{$node}{'RAxML_Number_Deleted'} = $no_deletions;
	#This is the number of traits created/deleted along the ancestral branch ABOVE the current node 
	
	unless ($TreeCacheHash->{$node}{'is_Leaf'}){
			
		my @Children = @{$TreeCacheHash->{$node}{'each_Descendent'}};
			
		$TreeCacheHash->{$node}{'RAxML_Total_Number_Created'} = $no_creations;
		$TreeCacheHash->{$node}{'RAxML_Total_Number_Deleted'} = $no_deletions;	
			
		foreach my $Child (@Children){
			
			 my ($Total_no_creations,$Total_no_deletions) = RAxML_Ancestral_Trait_Changes_in_Clade($TreeCacheHash,$Child,$root,$traitlabelsarray);
			#Notice the replacement of $node with $Child.
			#Total_no_X is the number of deletions/ creations in branches BELOW this node. 
			#This will be added to the number of deletions along this nodes ancestral branch to give the total in the clade BELOW the branch beneath the ancestor of the node under study

			$TreeCacheHash->{$node}{'RAxML_Total_Number_Created'} += $Total_no_creations;
			$TreeCacheHash->{$node}{'RAxML_Total_Number_Deleted'} += $Total_no_deletions;
		}
		
		my ($Total_no_creations,$Total_no_deletions) = ($TreeCacheHash->{$node}{'RAxML_Total_Number_Created'},$TreeCacheHash->{$node}{'RAxML_Total_Number_Created'});
		
		return($Total_no_creations,$Total_no_deletions);
			
	}else{
			
		$TreeCacheHash->{$node}{'RAxML_Total_Number_Created'} = $no_creations;
		$TreeCacheHash->{$node}{'RAxML_Total_Number_Deleted'} = $no_deletions;	
			
		return($no_creations,$no_deletions);	
	}
}


1;

