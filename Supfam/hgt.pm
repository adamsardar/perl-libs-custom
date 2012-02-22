package Supfam::hgt;
require Exporter;

=head1 NAME

Supfam::hgt

=head1 SYNOPSIS

Holds functions related to calculating infromation regarding horizontal gene transfer 
use Supfam::hgt;

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
				RandomModel
				RandomModelJulian
				DeletedJulian
				DeletedPoisson
				RandomModelPoisson
				RandomModelCorrPoisson
				RandomModelCorrPoissonOptimised
				calculatePosteriorQuantile
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;
use Bio::TreeIO;

use Bio::Tree::TreeFunctionsI;
use Time::HiRes;
use POSIX qw(floor ceil);

use Math::Random qw(random_poisson random_uniform random_uniform_integer);
use Set::IntervalTree; #Used in optimised code for mapping deletion events onto the phylogenetic tree

use lib "$ENV{HOME}/bin/perl-libs-custom";
use Supfam::TreeFuncsNonBP;
use Supfam::Utils;

use Data::Dumper;

sub DeletedJulian($$$$$$$){ 
	#($subtree,0,0,$NodesObserved)
	my ($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$treeroot,$DomArch) = @_;
	#$time in this function is evolutionary time and $dels are the number of observed deletions in that time
		
    my @InternalNodes = @{$TreeCacheHash->{$MRCA}{'each_Descendent'}};
    @InternalNodes = ($MRCA) if($TreeCacheHash->{$MRCA}{'is_Leaf'});
    #InternalNodes will be a list of all nodes within the tree, grown in the loop below. Some initial values are thrown in here.
    my @GenomesWithTraits = keys(%{$HashOfGenomesObserved});
      
    #  print Dumper($HashOfGenomesObserved);
      
       map{if($HashOfGenomesObserved->{$DomArch}{$_} > 1){print $DomArch."\n"; print $_."\n"; die '$HashOfGenomesObserved should have UNIQUE key value pairs. '.$_.' is non unique. There is a serious issue with the way that data is inputed into the hash'."\n";}}(keys(%{$HashOfGenomesObserved->{$DomArch}}));
       # Error check. Note that this sub requires that the two lists of genomes be unique (no internal repeats in the lists)!
		
	   while(scalar(@InternalNodes)){
	    	
	    	my $node = pop(@InternalNodes);
	    	
	    	my $DistanceFromAncestor = $TreeCacheHash->{$node}{'branch_length'};
 	
	    	my @NodeDescendents = map{$TreeCacheHash->{$_}{'node_id'}}@{$TreeCacheHash->{$node}{'Clade_Leaves'}}  ;# All descendent leaf ids
			@NodeDescendents = ($TreeCacheHash->{$node}{'node_id'}) if($TreeCacheHash->{$node}{'is_Leaf'}); # This allows for when the node is a leaf and hence has no descendents
	  	     	   	  	   
			my (undef,$Intersection,undef,undef) = IntUnDiff(\@NodeDescendents,\@GenomesWithTraits);		
					
			if(scalar(@$Intersection)){ # if there is any overlap between the list of genomes containing dom arch and those in the clade being studied
		
	    		$time += $DistanceFromAncestor;
	    		
	    		unless ($TreeCacheHash->{$node}{'is_Leaf'}){
	    			
					my @NodeChildren = @{$TreeCacheHash->{$node}{'each_Descendent'}};
	    			push(@InternalNodes,@NodeChildren);
	    		}
	 
	    	}else{
	    		
	    		$time += $DistanceFromAncestor/2;
	    		$dels++;
	    	}
	    }
	   
		return($dels,$time);
}

=pod
=item * DeletedJulian
This function looks at deletions of a particular domain archtecture of interest. The algorithm is as per Julian Gough's original hgt script. Given a tree and a list of genomes
 that posses that domain architecture, this function moves through the tree, level by level, clade by clade. If a clade contains an example of the domain arch, add the evolutionary
 time distance of the clade root from its parent and look at the clades of its children nodes. If no examples of the dom arch are found, increase the number of deletions by one and move on. 
Deleted($BioPerltreeobject,$CladeParentNode,$dels,$time,$GenomesPossessingDomainArchitecture)
=cut

sub RandomModelPoisson($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	
	my @PoissonianDeletions = random_poisson($Iterations,$Expected_deletions); #Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my $UniformDeletions = [];
	
	foreach my $DeletionSimultation (@PoissonianDeletions){ #For $Iterations
		
		@$UniformDeletions = random_uniform($DeletionSimultation,0,1); # Number of deletions, drawn from a poissonian above, uniformly distributed across the tree.
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		#$DeletionsNumberDistribution->{$DeletionSimultation}++;
		
		foreach my $DeletionPoint (@$UniformDeletions) {
                    
			my $index = 0; 
			while ($DeletionPoint > $ProbabilityIntervals[$index]){$index++;}
			# @ProbabilityIntervals is a precalculated hash of all the nodes in the sub-tree from the MRCA ($root) and where they sit in a stretched out sum of all branch lengths. $DeletedNode is a uniform random point along this line.
			#The above while loop is used to find the suitable point at which 
			my $DeletedNode = $ProbabilityHash->{$ProbabilityIntervals[$index]};
			
			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedNode}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedNode}) if ($TreeCacheHash->{$DeletedNode}{'is_Leaf'});		
		}
		
		my $no_model_genomes = scalar(keys(%ModelCladeGenomesHash));
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my $SelftestValue = $$RawResults[scalar(rand(@$RawResults))]; # A single uniform random simulation value
	
	#my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	#my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value
	
	return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution);
}

=pod
=item * RandomModelPoisson
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.
=cut

sub RandomModelCorrPoisson($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	
	my @PoissonianDeletions = random_poisson($Iterations,$Expected_deletions); #Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my $UniformDeletions = [];
	
	while  (@PoissonianDeletions){ #For $Iterations number of times
	
		my $DeletionSimultation = shift(@PoissonianDeletions); #remove the first entry in array and set it as the number of deletions in this simulation
		
		@$UniformDeletions = random_uniform($DeletionSimultation,0,1); # Number of deletions ($DeletionSimultation), drawn from a poissonian above, uniformly distributed across the tree.
		#Mallocing constantly
		
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		
		#$DeletionsNumberDistribution->{$DeletionSimultation}++;
		
		foreach my $DeletionPoint (@$UniformDeletions) {
                    
			my $index = 0; 
			while ($DeletionPoint > $ProbabilityIntervals[$index]){$index++;}
			# @ProbabilityIntervals is a precalculated hash of all the nodes in the sub-tree from the MRCA ($root) and where they sit in a stretched out sum of all branch lengths. $DeletedNode is a uniform random point along this line.
			#The above while loop is used to find the suitable point at which a deletion occurs
			my $DeletedNode = $ProbabilityHash->{$ProbabilityIntervals[$index]};
			
			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedNode}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedNode}) if ($TreeCacheHash->{$DeletedNode}{'is_Leaf'});		
		}
		
		## test to see if the simulation has resulted in an entire clade possessing a domain architecture and nothing else (i.e. which would make us find a new MRCA and a deletion rate of 0)
		#OR that the simulation has ended with no domain archtectures present anywhere OR that they are present everywhere, in which case we would set the deletion rate as zero
	
		my @ModelRemianingLeaves = keys(%ModelCladeGenomesHash);
		my $ModelFullCladeExclusive = 0; #Preallocate
		
		my $no_model_genomes = scalar(@ModelRemianingLeaves);
		
		if($no_model_genomes > 0){

			my $ModelRoot = FindMRCA($TreeCacheHash,$root,\@ModelRemianingLeaves);
			my @ModelFullCladeLeaves = @{$TreeCacheHash->{$ModelRoot}{'Clade_Leaves'}};
						
			(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff(\@ModelFullCladeLeaves,\@ModelRemianingLeaves)	; #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
		}
		
		if ($no_model_genomes == 0  || $no_model_genomes == scalar(@CladeGenomes) || scalar(@$ModelFullCladeExclusive) == 0){

			push(@PoissonianDeletions,random_poisson(1,$Expected_deletions)); #push a number onto the end on the deletions array
			next;
		} #IFF the simulation has ended with no genomes possesing the architecture (extinction) or with complete ubiquity in the clade under study,
		# or we have ubiquity in the clade beneath the MRCA of the simulated genomes
		# we discard the result (these three conditions would mean that we wouldn't be studying the domain architecture, leading to bias)
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value
	return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution);
}

=pod
=item * RandomModelCorrPoisson
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.

Unlike RandomModelPoisson, however, this implementation ignores results that would produce 0 (i.e. and extinct architecture) or
100% abndance (i.e. with a deletion rate of 0 and exluded from our models).
=cut


sub RandomModel {
	
	my ($CladeGenomes,$ObservedGenomes,$TreeCacheHash,$DeletionRate,$MRCA) = @_;
	
	#Note, these $*Genome variables are BioPerl nodeID objects
	
	my $LHHash={};
	$LHHash->{$MRCA}=1;
	
	my $TreeLikelihood = 1;

	my @UnobservedIDs;

	foreach my $CladeGenomeID (@$CladeGenomes){
		
		if(grep{/$CladeGenomeID/}@$ObservedGenomes){ #If $CladeGenome is in the list of genomes containing architecture
			
		my $ObservedGenomeLH = PositiveDomArchObservations($LHHash,$CladeGenomeID,$TreeCacheHash,$DeletionRate); #Populate $LHHash with LH values resulting from observed values
		$TreeLikelihood = $TreeLikelihood*$ObservedGenomeLH;
			
		}else{
			
			push(@UnobservedIDs,$CladeGenomeID);
			}
	}

	#print STDERR "Unobserved genome likelihoods:\n";
	
	foreach my $UnobservedGenome  (@UnobservedIDs){

		my $UnobservedGenomeLH = NegativeDomArchObservations($LHHash,$UnobservedGenome,$TreeCacheHash,$DeletionRate); #Populate $LHHash with LH values
		$TreeLikelihood = $TreeLikelihood*(1-$UnobservedGenomeLH); #The value is 1-prob as we are dealing with the probability of NOT observing dom arch

	}
	
	return($TreeLikelihood);

}

=pod
=head2 Methods
=over 4
=cut

sub RandomModelJulian($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;

	#$root is the root of the subtree or the most recent common ancestor
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	
    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength;
	
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my $UniformDeletions = [];
	
	for (1 .. $Iterations){ #For $Iterations
		
		my $DeletionSimultation = POSIX::floor($Expected_deletions);
		$DeletionSimultation+=1 if(rand() < ($Expected_deletions-$DeletionSimultation)); #Establises the number of deletions to be used in the simulation
		
		#$DeletionsNumberDistribution->{$DeletionSimultation}++;
		
		@$UniformDeletions = random_uniform($DeletionSimultation); # Draw uniform(0,1) numbers so that we have an array of length equal to the number of deletions expected
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		
		foreach my $DeletionPoint (@$UniformDeletions) {
                    
			my $index = 0; while ($DeletionPoint > $ProbabilityIntervals[$index]){$index++;}
			# @ProbabilityIntervals is a precalculated hash of all the nodes in the sub-tree from the MRCA ($root) and where they sit in a stretched out sum of all branch lengths. $DeletedNode is a uniform random point along this line.
			my $DeletedNode = $ProbabilityHash->{$ProbabilityIntervals[$index]};
			
			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedNode}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedNode}) if ($TreeCacheHash->{$DeletedNode}{'is_Leaf'});		
		}
		
		my $no_model_genomes = scalar(keys(%ModelCladeGenomesHash));
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value
	
	return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution);
	
}


sub RandomModelJulianOpt($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;

	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength;
	
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	for (1 .. $Iterations){ #For $Iterations
		
		my $DeletionSimultation = POSIX::floor($Expected_deletions);
		$DeletionSimultation+=1 if(rand() < ($Expected_deletions-$DeletionSimultation));
		
		my @UniformDeletions = random_uniform($DeletionSimultation); # Number of deletions, drawn from a poissonian above, uniformly distributed across the tree.
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		
		foreach my $DeletionPoint (@UniformDeletions) {
                    
			my $index = 0; while ($DeletionPoint > $ProbabilityIntervals[$index]){$index++;}
			# @ProbabilityIntervals is a precalculated hash of all the nodes in the sub-tree from the MRCA ($root) and where they sit in a stretched out sum of all branch lengths. $DeletedNode is a uniform random point along this line.
			my $DeletedNode = $ProbabilityHash->{$ProbabilityIntervals[$index]};
			
			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedNode}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedNode}) if ($TreeCacheHash->{$DeletedNode}{'is_Leaf'});		
		}
		
		my $no_model_genomes = scalar(keys(%ModelCladeGenomesHash));
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value

	return($SelftestValue,$distribution,$RawResults);
}


=pod
=item * RandomModelPoisson
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.
=cut

sub calculatePosteriorQuantile($$$$){

	my ($SingleValue,$DistributionHash,$NumberOfSimulations,$CladeSize) = @_;
		
	$DistributionHash->{$SingleValue}++;
	my $NumberOfSimulationsLT = 0;
	
	my @DistributionIndicies = keys(%$DistributionHash);

	foreach my $Number_of_genomes  (0 .. $CladeSize){
		
		last if($Number_of_genomes == $SingleValue);
		
		if (exists($DistributionHash->{$Number_of_genomes})){
			
			$NumberOfSimulationsLT += $DistributionHash->{$Number_of_genomes};
		}
   }
	
	my $Degeneracy = $DistributionHash->{$SingleValue};#Number of simulations of equal score. We place our point to sum up to uniform in this region
	my ($DegeneracyContribution) = random_uniform(1,0,$Degeneracy);
	$NumberOfSimulationsLT += $DegeneracyContribution;
	#$NumberOfSimulationsLT += $Degeneracy/2;
	
	my $PosteriorQuantile = $NumberOfSimulationsLT/$NumberOfSimulations;
	
	$DistributionHash->{$SingleValue}--;
	#Undo the modification to the distribution
	
	return($PosteriorQuantile);
}


=pod * calculatePosteriorQuantile($SingleValue,$DistributionHash)

Calculates where in the total area in a probability distribution a single value occurs. Returns a value between 0 and 1. These should be uniformly distributed, from the simple fact that sum(andy distribution)
=1 and we are choosing a unifrom point from theis area.
=cut


sub PositiveDomArchObservations{
	#Removed the prototyping - be careful

	
	my ($LHHash,$NodeID,$TreeCacheHash,$DeletionRate) = @_;
	my $NodeLH;
	my $ParentNodeID = $NodeID->ancestor;
	my $ParentLH; #This is the probability that the domain architecutre under study WAS present in the parent genome, given previous evidence.
	
	unless(exists($LHHash->{$ParentNodeID})){
	
		$ParentLH = PositiveDomArchObservations($LHHash,$ParentNodeID,$TreeCacheHash,$DeletionRate);
		$LHHash->{$ParentNodeID} = 1; # This is a consequence of dollo parsimony - the fact that there was already a presence in the LHHash means that the probability of the ancestral state has already been factored in
		
	}else{
		
		$ParentLH = $LHHash->{$ParentNodeID}; #This will always be one, see above comment. Note that this is not the case in NegativeDomArchObservations sub
	}
		
	my $BranchLength = $TreeCacheHash->{$NodeID}{'branch_length'};
	
	my $ProbOfDeletion = 1-exp(-$DeletionRate*$BranchLength); # The probability of a deletion on the branch between ancestor and given node at a given deltion rate. Based on the cumulative density function of the exponetial distribution
	
	$NodeLH = $ParentLH*(1 - $ProbOfDeletion); #i.e. probability of DA being present in ancestor given observations multiplied by the probability of dom arch NOT being deleted
	$LHHash->{$NodeID} = $NodeLH; #Update hash with the result. Note that this is the probability of observing the domain arch in the given node given the previously observed nodes. This is NOT an independant likelihood, but a marginal
	
	return($NodeLH);
}


=pod * PositiveDomArchObservations
This is a sub to calculate the likelihood of a node containing the domain architecture in question, given previous observations (stored in $LHHash).
=cut

sub NegativeDomArchObservations{
		#Removed the prototyping - be careful
	
	my ($LHHash,$NodeID,$TreeCacheHash,$DeletionRate) = @_;
	my $NodeLH;
	my $ParentNodeID = $NodeID->ancestor;
	my $ParentLH; #This is the probability that the domain architecutre under study WAS present in the parent genome, given previous evidence.
	
	unless(exists($LHHash->{$ParentNodeID})){
	
		$ParentLH = PositiveDomArchObservations($LHHash,$ParentNodeID,$TreeCacheHash,$DeletionRate);
					
	}else{
		
		$ParentLH = $LHHash->{$ParentNodeID};
	}
		
	my $BranchLength = $TreeCacheHash->{$NodeID}{'branch_length'};
	
	my $ProbOfDeletion = 1-exp(-$DeletionRate*$BranchLength); # The probability of a deletion on the branch between ancestor and given node at a given deltion rate. Based on the cumulative density function of the exponetial distribution
	
	$NodeLH = $ParentLH*(1 - $ProbOfDeletion); #i.e. probability of DA being present in ancestor given observations multiplied by the probability of dom arch NOT being deleted
	$LHHash->{$NodeID} = $NodeLH; #Update hash with the result. Note that this is the probability of observing the domain arch in the given node given the previously observed nodes. This is NOT an independant likelihood, but a marginal
	
	return($NodeLH);
}

=pod * NegativeDomArchObservations
This is a sub to calculate the likelihood of a node containing the domain architecture in question, given previous observations (stored in $LHHash). Note the difference between these two functions:
One assumes that the inputed node has had the domain architecture observed in its genomes. Using strict dollo parsimony therefore, all later references to whther an ancestor contains that dom arch will
be one, as it is the only possible way to explain the previously observed sighiting of the domain architecture.
=cut

1;
