package Supfam::hgt;
require Exporter;

=head1 NAME

Supfam::hgt

=head1 SYNOPSIS

Holds functions related to calculating and processing information regarding horizontal gene transfer of domain architectures
use Supfam::hgt;

=head1 AUTHOR

Adam Sardar (adam.sardar@bris.ac.uk)

=head1 COPYRIGHT

Copyright 2011-2012 Gough Group, University of Bristol.

=head1 SEE ALSO

Supfam::Config.pm

=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
				DeletedJulian
				DeletedJulianDetailed
				HGTTreeDeletionModel
				HGTTreeDeletionModelOptimised
				RandomModelCorrPoissonDeletionDetailed
				RandomModelCorrPoissonOptimised
				RandomModelCorrPoissonOptimisedDetailed
				calculatePosteriorQuantile
				calculateOldStylePosteriorQuantile
				calculateJulianStylePosteriorQuantile
				calculateHashContinuousPosteriorQuantile
				calculateContinuousPosteriorQuantile
				RandomModelPoissonOptimised
				HGTshuffle
				DeletedSimAnneal
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use lib "$ENV{HOME}/bin/perl-libs-custom/";

use Modern::Perl;

use Time::HiRes;
use POSIX qw(floor ceil);

use Math::Random qw(random_poisson random_uniform random_uniform_integer random_exponential random_negative_binomial random_gamma);

use Supfam::TreeFuncsNonBP;
use Supfam::Utils;
use Supfam::PointTree; #For use in the optimised CorrPoiss

use Data::Dumper;
use List::Util qw(sum);
use Math::Decimal qw(dec_cmp is_dec_number);
use Algorithm::Combinatorics qw(combinations);
use Carp;
use Carp::Assert;
use Carp::Assert::More;
use List::Compare;

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

sub DeletedJulianDetailed($$$$$$){ 
	
	my ($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$treeroot,$DomArch) = @_;
	#$time in this function is evolutionary time and $dels are the number of observed deletions in that time
	
	my $DeletionPoints = []; #This will be a list of nodes that are direct decendents to deletion points
	
    my @InternalNodes = @{$TreeCacheHash->{$MRCA}{'each_Descendent'}};
    @InternalNodes = ($MRCA) if($TreeCacheHash->{$MRCA}{'is_Leaf'});
    #InternalNodes will be a list of all nodes within the tree, grown in the loop below. Some initial values are thrown in here.
    my @GenomesWithTraits = keys(%{$HashOfGenomesObserved});
      
       map{if($HashOfGenomesObserved->{$_} > 1){print $DomArch."\n"; print $_."\n"; die '$HashOfGenomesObserved should have UNIQUE key value pairs. '.$_.' is non unique. There is a serious issue with the way that data is inputed into the hash'."\n";}}(keys(%{$HashOfGenomesObserved}));
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
	    		
	    		push(@$DeletionPoints,$node);
	    	}
	    }
	   
	   #calculate the distances between deletion points now.
	   
	   #Find MRCA of two points
	   #Take two leaves, one from the clades beneath each deletion
	  
	   #Calculate the distance between these two points
	   #Write a sub that sums ancestral distance to MRCA. Then add half the distance above branch
	   
	   my $InterDeletionDistances = [];
	   
	   
	   unless(scalar(@$DeletionPoints) < 2){
	   	
			my $iter = combinations(\@$DeletionPoints,2);
			#Caculate all n choose 2 combinations. We shall study the distances between each point
			
			while (my $combination = $iter->next) {
			   
			  my ($DelPointA,$DelPointB) = @$combination;
			  
			  my $ArrayRefA = $TreeCacheHash->{$DelPointA}{'Clade_Leaves'};
			  $ArrayRefA = [$DelPointA] if($TreeCacheHash->{$DelPointA}{'is_Leaf'} == 1);
			  
			  my $ArrayRefB = $TreeCacheHash->{$DelPointB}{'Clade_Leaves'};
			  $ArrayRefB = [$DelPointB] if($TreeCacheHash->{$DelPointB}{'is_Leaf'} == 1);  
			  
			  my $UnionHash = {};
			  map{$UnionHash->{$_}=undef}@$ArrayRefA;
			  map{$UnionHash->{$_}=undef}@$ArrayRefB;
			  #Use a hash to quickly find the union of all the genomes in the two lists.
			  
			  my @Union = keys(%$UnionHash);
			   
			  my $MRCA = FindMRCA($TreeCacheHash,$treeroot,\@Union);
			  
			 # print $MRCA."\n";
			 # print " => MRCA\n";
			  
			  my $InterDeletionPoint = 0;
			  
			  my $ABranchSum = MRCADistanceSum($TreeCacheHash,$DelPointA,$MRCA);
			  $InterDeletionPoint += $ABranchSum;
			  
			  my $BBranchSum = MRCADistanceSum($TreeCacheHash,$DelPointB,$MRCA);
			  $InterDeletionPoint += $BBranchSum;
			  #Add the distance from the nodes from MRCA to the sum
			   
			  $InterDeletionPoint -= ($TreeCacheHash->{$DelPointA}{'branch_length'}/2);
			  $InterDeletionPoint -= ($TreeCacheHash->{$DelPointB}{'branch_length'}/2);
				#subtract half the branch length between nodes and their direct ancestor. This is because we don't know exactly where on the branch this occured.
			
			  push(@$InterDeletionDistances,$InterDeletionPoint);
			  #$InterDeletionPoint is the distance between two deletion points. We assume that a deletion occurs midway on a branch, hence the subtraction of the values from the total sum of the points distance from the MRCA
			}
	   }

	   #Push onto array of distances
	   return($dels,$time,$InterDeletionDistances);
}

=pod
=item * DeletedJulianDetailed($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$treeroot,$DomArch)

TThis function looks at deletions of a particular domain archtecture of interest. The algorithm is as per Julian Gough's original hgt script. Given a tree and a list of genomes
 that posses that domain architecture, this function moves through the tree, level by level, clade by clade. If a clade contains an example of the domain arch, add the evolutionary
 time distance of the clade root from its parent and look at the clades of its children nodes. If no examples of the dom arch are found, increase the number of deletions by one and move on. 

this is a more vebose version of  DeletedJulian. It calculates where all deletions have occured before calculating all n choose 2 distances between them.

=cut


sub DeletedSimAnneal{ 

	my ($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$treeroot,$DomArch,$CladeSize,$opt_iterations,$HGTpercentage,$model,$NoGenomesObserved) = @_;
	#$time in this function is evolutionary time and $dels are the number of observed deletions in that time
	

	#s ← s0; e ← E(s)                                  // Initial state, energy.
	#sbest ← s; ebest ← e                              // Initial "best" solution
	#k ← 0                                             // Energy evaluation count.

	my $SimIterations = 100;
	
	my $lambda_original = $dels/$time;
	my $lamba_best = $lambda_original;
	
	my (undef,$predistribution,undef,undef) = HGTTreeDeletionModelOptimised($MRCA,$model,$SimIterations,[$lamba_best],$TreeCacheHash,$HGTpercentage/100);
	
	$predistribution->{$NoGenomesObserved}++; #Prevents perl from kicking out by a divide by zero error here. A pseudocount
	
	my $LH = $predistribution->{$NoGenomesObserved};
	my $BestLambda_Energy = $LH/($SimIterations+1);

	my $OriginalEnergy = $BestLambda_Energy;
	
	#print $BestLambda_Energy;
	
	my $count = 1;
	my ($gamma_alpha,$gamma_beta) = ($dels,$time);

		
	while($count < $opt_iterations){
		
	
		my $Temperature = $opt_iterations/$count -1; #This approaches 0, steeply at first and then slows out. Diff is 1/count^2
		
		my $lambda_new = random_gamma(1,$gamma_beta,$gamma_alpha);
		#Choose a new value for lambda from the vaccinity of the current best values

		my (undef,$simdistribution,undef,undef) = HGTTreeDeletionModelOptimised($MRCA,$model,$SimIterations,[$lambda_new],$TreeCacheHash,0);
		$simdistribution->{$NoGenomesObserved}++; #Prevents perl from kicking out by a divide by zero error here. A pseudocount
		my $LH = $simdistribution->{$NoGenomesObserved};
		my $NewLambda_Energy = $LH/($SimIterations+1);
		#Self test treats a randomly chosen simulation as though it were a true result. We therefore reduce the distribution count at that point by one, as we are picking it out. This is a sanity check.
					
		if(random_uniform() <= exp(($NewLambda_Energy - $BestLambda_Energy)/$Temperature)){
			
			#Move vacinity - perform a single simulation at our new lambda rate and then measure the number of deletions.
			my $SingleCombGenomeSimHash = HGTTreeDeletionModelOptimised($MRCA,$model,1,[$lambda_new],$TreeCacheHash,0);
			my $SingleSimDomCombGenomeHash = {};
			map{$SingleSimDomCombGenomeHash->{'Comb'}{$TreeCacheHash->{$_}{'node_id'}}=1}(keys(%$SingleCombGenomeSimHash));
			
			my @SimGens = keys(%{$SingleCombGenomeSimHash});
			
			#print join("\t",@SimGens);
			#print "\n";
			
			my $submrca;
			$submrca = $SimGens[0];
			
			unless(scalar(@SimGens) <= 1){
			
				$submrca = FindMRCA($TreeCacheHash,$treeroot,\@SimGens) ;#($TreeCacheHash,$root,$LeavesArrayRef);
				
				($gamma_alpha, $gamma_beta) = DeletedJulian($submrca,0,0,$SingleSimDomCombGenomeHash,$TreeCacheHash,$treeroot,'Comb'); # ($tree,$AncestorNodeID,$dels,$time,$GenomesOfDomArch) - calculate deltion rate over tree
				#Update the 'vaccinity paramters'. These are the paramters that feed into the gamma rate choice step	
			}
			
		}
		
		
		if($NewLambda_Energy > $BestLambda_Energy){
			
			$lamba_best = $lambda_new;
			$BestLambda_Energy = $NewLambda_Energy;
		}
		
		$count++;
	}

	my ($selftest,$Distribution,undef,undef) = HGTTreeDeletionModelOptimised($MRCA,$model,$SimIterations,[$lamba_best],$TreeCacheHash,0);	
	
	$Distribution->{$selftest}-- if(exists($Distribution->{$selftest}));
	my $selftestval = calculatePosteriorQuantile($NoGenomesObserved,$Distribution,$SimIterations+1,$CladeSize);
	
	return($lamba_best, $lambda_original,$BestLambda_Energy,$OriginalEnergy,$selftestval);

}

=pod
=item * DeletedJulian

A new and crazy idea of mine - use simmualted annealing to maximise the likelihood of the parameter lambda (rate parameter) in the poisson model. I have literally no
idea how long this will take or even if this is going to be  numerically fast enough. Lets find out ....

=cut



sub HGTTreeDeletionModel($$$$$$$) {
	
	my ($root,$model,$Iterations,$ndelsobs,$timedelsobserved,$TreeCacheHash,$HGTpercentage) = @_;
	#$root is the root of the subtree or the most recent common ancestor

	my $corrflag = ($model =~ m/corr/i)?1:0;
	
    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash; map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	
	my $deletion_rate = $ndelsobs/$timedelsobserved;
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	#This is the average number of deletions expected at a given deletion rate. This is the single paramentr of input into a poisson model
	
	####
	
	my @NumberOfDeletions;
	
	if($model =~ m/poisson/i){
		
		@NumberOfDeletions = random_poisson($Iterations,$Expected_deletions);
		#Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE for the exponential formulation of the poisson process)
				
	}elsif($model =~ m/negbin/i){
		
		@NumberOfDeletions = random_negative_binomial($Iterations, $ndelsobs, $TotalBranchLength/($TotalBranchLength + $timedelsobserved));
		#Number of deletions in this iteration. This is drawn from a negative binomial distribution with parameters determined so as to ensure that the formulation as a gamma-poisson mixture continues.
	}
	#####
	
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my $UniformDeletions = [];
	
	my $SingleSimGenomeHash = {};
	my @HGTUniformSimsPool = random_uniform($Iterations,0,1);
	#Crete a pool of unifrom random numbers for determining if an HGT has occured or not.
	
	while (@NumberOfDeletions){
		
		my $DeletionSimultation = pop(@NumberOfDeletions);
		
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
		
		if(pop(@HGTUniformSimsPool) < $HGTpercentage){
			#Random assignment of dom arches to genomes
			
			my @GenomesWithDAByHGT = @{HGTshuffle(\@CladeGenomes,'Power')};
			#At current we only allow for the current model to be used.
			
			if(scalar(@GenomesWithDAByHGT) > 0){
				map{$ModelCladeGenomesHash{$_} =1;}@GenomesWithDAByHGT;
				#Add a chunk of domain architectures into the clade genomes list
			}
			
		}
			
		my @ModelRemianingLeaves = keys(%ModelCladeGenomesHash);
		my $ModelFullCladeExclusive = 0;
		
		my $no_model_genomes = scalar(@ModelRemianingLeaves);
		
		if($corrflag){
			
			if($no_model_genomes > 0){
	
				my $ModelRoot = FindMRCA($TreeCacheHash,$root,\@ModelRemianingLeaves);
				my @ModelFullCladeLeaves = @{$TreeCacheHash->{$ModelRoot}{'Clade_Leaves'}};
							
				(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff(\@ModelFullCladeLeaves,\@ModelRemianingLeaves)	; #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
			}
			
			if ($no_model_genomes == 0  || $no_model_genomes == scalar(@CladeGenomes) || scalar(@$ModelFullCladeExclusive) == 0){
	
				if($model =~ m/poisson/i){
					
					@NumberOfDeletions = random_poisson($Iterations,$Expected_deletions);		
				}elsif($model =~ m/negbin/i){
					
					@NumberOfDeletions = random_negative_binomial($Iterations, $ndelsobs, $TotalBranchLength/($TotalBranchLength + $timedelsobserved));
		
				}#push a number onto the end on the deletions array
				
				push(@HGTUniformSimsPool,random_uniform(1,0,1)); #push a number onto the end on the uniform pool array
				next;
			} 
			
			#IFF the simulation has ended with no genomes possesing the architecture (extinction) or with complete ubiquity in the clade under study,
			# or we have ubiquity in the clade beneath the MRCA of the simulated genomes
			# we discard the result (these three conditions would mean that we wouldn't be studying the domain architecture, leading to bias)
		}
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	
		if ($Iterations == 1){
			
			%$SingleSimGenomeHash = %ModelCladeGenomesHash;
		}

	}
	
	my $SelftestValue = $$RawResults[scalar(rand(@$RawResults))]; # A single uniform random simulation value
	
	#my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	#my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value
	
	unless($Iterations == 1){
		
		return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution);
		
	}else{
		
		return($SingleSimGenomeHash);
	#Allows for a single simulation of the model to be performed and dumped out
	}
}

=pod
=item * HGTTreeDeletionModel ($root,$model,$Iterations,$ndelsobs,$timedelsobserved,$TreeCacheHash,$HGTpercentage)
A deletion model based on the poisson/negativ binomial distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.
=cut


sub HGTTreeDeletionModelOptimised {
	
	my ($root,$model,$Iterations,$modelparams,$TreeCacheHash,$HGTprob,$HGTmodel) = @_;
	#$root is the root of the subtree or the most recent common ancestor, HGTmodel is the type of model 
	
	assert_in($model,['corrpoisson','corrnegbin','poisson','negbin']);
	assert_in($HGTmodel,['scatter','drop'],'Allowed HGTmodels are "scatter" and "drop"');
	#Ensure that the models are set up correctly
	
	my ($deletion_rate,$ndelsobs,$timedelsobserved);
	
	if (scalar(@$modelparams) == 2){
		
		($ndelsobs,$timedelsobserved) = @$modelparams;
		$deletion_rate = $ndelsobs/$timedelsobserved;
		
	}elsif(scalar(@$modelparams) == 1){
		
		($deletion_rate) = @$modelparams;
		
	}else{
		
		carp "Too many model paramters provided. WTF!?\n";
	}
	
	assert(scalar(@_) == 7, "This function takes 7 arguments\n");

	my $corrflag = ($model =~ m/corr/i)?1:0;
	
    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash; map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree

	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	#This is the average number of deletions expected at a given deletion rate. This is the single paramentr of input into a poisson model
	
	##

    my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
    my @Points = keys(%$ProbabilityHash);
        
    my $PointTree = Supfam::PointTree->new;
    $PointTree->build(\@Points);
  	#Create the PointTree - this is for fast searching of the linear space to find where a deletion ahs occured on our tree
   
	my @NumberOfDeletions;
	
	if($model =~ m/poisson/i){
		
		@NumberOfDeletions = random_poisson($Iterations,$Expected_deletions);
		#Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE for the exponential formulation of the poisson process)
				
	}elsif($model =~ m/negbin/i){
		
		carp "Issue here - ndels is undefined!" if ($ndelsobs ~~ undef);
		@NumberOfDeletions = random_negative_binomial($Iterations, $ndelsobs, $TotalBranchLength/($TotalBranchLength + $timedelsobserved));
		#Number of deletions in this iteration. This is drawn from a negative binomial distribution with parameters determined so as to ensure that the formulation as a gamma-poisson mixture continues.
	}
	#####
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	my $DetailedRawResults = [];
	my $SingleSimGenomeHash = {};
	#Prealloacte a few variables used to hold simulation results
	
	my @HGTUniformSimsPool = random_uniform(2*$Iterations,0,1) if($HGTprob > 0);
	#Crete a pool of unifrom random numbers for determining if an HGT has occured or not.
	
	my $NumberOfDelPoints = ceil(0.2*$Iterations*$Expected_deletions);

	$PointTree->UniformAssign($NumberOfDelPoints);
	$PointTree->UniformPoolAssign($NumberOfDelPoints);
	## Uniformly set deletion points across the subtree of interest.
	## Also, create a pool of uniform random numbers
	
	while (@NumberOfDeletions){
		
		my $DeletionSimultation = pop(@NumberOfDeletions);
		
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		my $DeletionPoints = [];
		$PointTree->UniformDraw($DeletionSimultation,$DeletionPoints);
		# A number ($DeletionSimultation) of uniform random deletions, uniform randomly scattered over the tree.
	
		while (my $DeletionPoint = pop(@$DeletionPoints)){
			
			my $DeletedBranch = $ProbabilityHash->{$DeletionPoint};

			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedBranch}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedBranch}) if ($TreeCacheHash->{$DeletedBranch}{'is_Leaf'});		
		}
		#Assign these deleitons to the 'Model Genomes Hash' (i.e. delete all the descendent possessions of the DA)


		if($HGTprob){
			#A time saving excercise - if HGTPercentage is off, we can jsut ignore all the code to do with HGT simulaton

			my $NumberOfGenomesWithHGTUndiscernable = 0;
			my $NumberOfGenomesWithHGTDiscernable = 0;
			#Sum from a random leaf to root of tree

			while($PointTree->UniformPoolDraw < $HGTprob){
			
				#Follow a binomial distribution of HGT probablities.
				
				if($HGTmodel eq 'drop'){
					
					my $HGTPoint = [];
					$PointTree->UniformDraw(1,$HGTPoint);
					my $HGTBranch = $ProbabilityHash->{$$HGTPoint[0]};
					
					my $HGTBranchToRootDistance = MRCADistanceSum($TreeCacheHash,$HGTBranch,$root);
					#Sum from HGT ancestor to root of tree
					
					unless($TreeCacheHash->{$HGTBranch}{'is_Leaf'}){
						
						my $lc = List::Compare->new('--unsorted', $TreeCacheHash->{$HGTBranch}{'Clade_Leaves'}, [keys(%ModelCladeGenomesHash)]);
						#Trying out a module from CPAN for computing the interesection of two lists
						
						my @GenomeIntersection = $lc->get_intersection;
						$NumberOfGenomesWithHGTUndiscernable += scalar(@GenomeIntersection);
						#Genome intersection that tells us if an HGT has resulted in an addition of an architecture to a point on the tree that already had the architecture or a 'new' part
						
						unless($NumberOfGenomesWithHGTUndiscernable){
												
							my %HGTSubTreeHash;
							map{$HGTSubTreeHash{$_} = 1; $ModelCladeGenomesHash{$_} = 1;}@{$TreeCacheHash->{$HGTBranch}{'Clade_Leaves'}};
							#Initialise an equivilent of 'Model Genome Hash', as with other simulations, but for this subtree only
							#Also, set all the leaves of the descendents of the HGT ancestor to be 1 in the full modelhash
																									
							my $SubTreeProbabilityHash = $TreeCacheHash->{$HGTBranch}{'Probability_Hash'};
		   					my @SubTreePoints = keys(%$SubTreeProbabilityHash);
		  					my $HGTSubPointTree = Supfam::PointTree->new;
		   					$HGTSubPointTree->build(\@SubTreePoints);
		   					#We're going to do a round of deletion simulations here, so construct a new set of objects for use in the simulation
		   					
							my $SubTreeTotalBranchLength = $TreeCacheHash->{$HGTBranch}{'Total_branch_lengths'};
							my $SubTreeExpected_deletions = $deletion_rate*$SubTreeTotalBranchLength;
							my ($SubTreeDeletions) = random_poisson(1,$SubTreeExpected_deletions);
							#A number of deletions, equal to the poisson distrbution over the sub tree at a given rate

							$DeletionPoints = [];
							$HGTSubPointTree->UniformDraw($SubTreeDeletions,$DeletionPoints);
							
							while (my $DeletionPoint = pop(@$DeletionPoints)){
			
								my $DeletedBranch = $SubTreeProbabilityHash->{$DeletionPoint};
								
								map{delete($ModelCladeGenomesHash{$_}) ; delete($HGTSubTreeHash{$_}) ;}@{$TreeCacheHash->{$DeletedBranch}{'Clade_Leaves'}};
								if ($TreeCacheHash->{$DeletedBranch}{'is_Leaf'}){delete($ModelCladeGenomesHash{$DeletedBranch}); delete($HGTSubTreeHash{$DeletedBranch}) ;}
								#For each deletion point, delete all descendents beneath that ancestor
							}
							#Run the deletion model on sub tree at the same rate as  before
							
							$NumberOfGenomesWithHGTDiscernable += scalar(keys(%HGTSubTreeHash));
						}
							
					}else{
						
						$ModelCladeGenomesHash{$HGTBranch} =1;
						#If the architecture has been transfered into a leaf node, set it as existing in %ModelCladeGenomesHash. Note that it could have already existed in 
					}
									
				}elsif($HGTmodel eq 'scatter'){
				
					#Random assignment of dom arches to genomes
					my @GenomesWithDAByHGT = @{HGTshuffle(\@CladeGenomes,'Power')};
					#At current we only allow for the current model to be used.
					
					if(scalar(@GenomesWithDAByHGT) > 0){
						map{$ModelCladeGenomesHash{$_} =1;}@GenomesWithDAByHGT;
						#Add a chunk of domain architectures into the clade genomes list
					}
				}
			}
		
			my $DetailedRawResultsString = join(":",(scalar(keys(%ModelCladeGenomesHash)),$NumberOfGenomesWithHGTUndiscernable,$NumberOfGenomesWithHGTDiscernable));
			push(@$DetailedRawResults,$DetailedRawResults);
		}
			
		my @ModelRemianingLeaves = keys(%ModelCladeGenomesHash);
		my $ModelFullCladeExclusive = 0;
		
		my $no_model_genomes = scalar(@ModelRemianingLeaves);
		
		if($corrflag){
			
			if($no_model_genomes > 0){
	
				my $ModelRoot = FindMRCA($TreeCacheHash,$root,\@ModelRemianingLeaves);
				my @ModelFullCladeLeaves = @{$TreeCacheHash->{$ModelRoot}{'Clade_Leaves'}};
							
				(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff(\@ModelFullCladeLeaves,\@ModelRemianingLeaves); #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
			}
			
			if ($no_model_genomes == 0  || $no_model_genomes == scalar(@CladeGenomes) || scalar(@$ModelFullCladeExclusive) == 0){
	
				if($model =~ m/poisson/i){
					
					@NumberOfDeletions = random_poisson(1,$Expected_deletions);		
				}elsif($model =~ m/negbin/i){
					
					@NumberOfDeletions = random_negative_binomial(1, $ndelsobs, $TotalBranchLength/($TotalBranchLength + $timedelsobserved));
		
				}#push a number onto the end on the deletions array
				
				pop(@$DetailedRawResults) if(@$DetailedRawResults); #Remove the last value from the HGT simualtion run, if indeed we're collecting the data
				next;
			} 
			
			#IFF the simulation has ended with no genomes possesing the architecture (extinction) or with complete ubiquity in the clade under study,
			# or we have ubiquity in the clade beneath the MRCA of the simulated genomes
			# we discard the result (these three conditions would mean that we wouldn't be studying the domain architecture, leading to bias)
		}
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		
		#Update the distribution of the run accordingly and store results in rawresults
	
		if ($Iterations == 1){
			
			%$SingleSimGenomeHash = %ModelCladeGenomesHash;
		}

	}
	
	my $SelftestValue = $$RawResults[scalar(rand(@$RawResults))]; # A single uniform random simulation value
	
	unless($Iterations == 1){
		
		return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution,$DetailedRawResults);
		
	}else{
		
		return($SingleSimGenomeHash);
	#Allows for a single simulation of the model to be performed and dumped out
	}
}

=pod
=item * HGTTreeDeletionModel ($root,$model,$Iterations,$ndelsobs,$timedelsobserved,$TreeCacheHash,$HGTpercentage)
A deletion model based on the poisson/negativ binomial distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.
=cut

sub HGTshuffle($$){
	
	my ($CladeGenomes,$delmodel) = @_;
										
			my $ShuffledCladeGenomes = [];
			@$ShuffledCladeGenomes = @{$CladeGenomes};
			
			my $RandomCladeInt;
			
			if($delmodel eq 'Uniform'){
				#Single unifrom number 'N' between 1 and size_of_clade.
				$RandomCladeInt = random_uniform_integer(1,0,scalar(@$ShuffledCladeGenomes)-1);
				
			}elsif($delmodel eq 'Power'){
				
				
				while($RandomCladeInt ~~ undef || $RandomCladeInt > (scalar(@$ShuffledCladeGenomes)-1) || $RandomCladeInt < 0){

					$RandomCladeInt = int(exp(-log(1-random_uniform(1,0,1))/1.25));
					
					#Alpha = 1.25 - paraemter of power law, results in an average of 5
					
					#Convert a uniform value to a power law distributed value by inverse transform sampling -> take a uniform random number and ask where it falls within the CDF of the power law
					#alpha = 1.25 was found through least squares regression with R> 0.95 in both eukaryotes and bacteria. p-value < 0.001
				}
				
				
			}elsif($delmodel eq 'Geometric'){
				
				#Single unifrom number 'N' between 1 and size_of_clade.
				while($RandomCladeInt ~~ undef || $RandomCladeInt > (scalar(@$ShuffledCladeGenomes)-1)){
					
					$RandomCladeInt = int(random_exponential(1,'12.29'));
					#12.29 is the average number of genomes that a dom arch belongs to in eukaryotes.
				}
			}
			
			#Shuffle the genomes, then choose the first 'N' terms.
			fisher_yates_shuffle($ShuffledCladeGenomes);
			
			my $HGTPossesions = [];
			push(@$HGTPossesions,@$ShuffledCladeGenomes[0 .. $RandomCladeInt]);
			#i.e. genomes that posses a DA as a consequence of HGT
			
	return($HGTPossesions);
}

=pod
=item * HGTshuffle
A very simple little model of HGT. A domain architecture can belong to one or more members of a list of genomes. This number is set to be equal to the number drawn from
a distribution of choice. N.B these distributions are parameterised by values drawn from Eukaryotes. Perhaps these will become options that you can set in future,

=cut

sub RandomModelPoissonOptimised($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
    
    my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
    my @Points = keys(%$ProbabilityHash);
        
    my $PointTree = Supfam::PointTree->new;
    $PointTree->build(\@Points);
  	#Create the PointTree
	
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	
	my @PoissonianDeletions = random_poisson($Iterations,$Expected_deletions); #Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	foreach my $DeletionSimultation (@PoissonianDeletions){ #For $Iterations
		
		my @UniformDeletions = random_uniform($DeletionSimultation,0,1); # Number of deletions, drawn from a poissonian above, uniformly distributed across the tree.
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		$DeletionsNumberDistribution->{$DeletionSimultation}++;
		
		foreach my $DeletionPoint (@UniformDeletions) {
                    
			my $index = 0; 
			my $DeletedNode = $PointTree->search($DeletionPoint);
			
			#print $DeletedNode."\n";
			
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


sub RandomModelCorrPoissonDeletionDetailed($$$$$) {
	
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
	my $UniformDeletions = [];
	
	my $SimulatedInterDeletionDistances = [];
	
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
		
		my $ModelRoot;
		
		if($no_model_genomes > 0){
			
			$ModelRoot = FindMRCA($TreeCacheHash,$root,\@ModelRemianingLeaves);
			
			my @ModelFullCladeLeaves = @{$TreeCacheHash->{$ModelRoot}{'Clade_Leaves'}};
						
			(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff(\@ModelFullCladeLeaves,\@ModelRemianingLeaves)	; #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
		}
		
		if ($no_model_genomes == 0  || $no_model_genomes == scalar(@CladeGenomes) || scalar(@$ModelFullCladeExclusive) == 0){

			push(@PoissonianDeletions,random_poisson(1,$Expected_deletions)); #push a number onto the end on the deletions array
			next;
		} #IFF the simulation has ended with no genomes possesing the architecture (extinction) or with complete ubiquity in the clade under study,
		# or we have ubiquity in the clade beneath the MRCA of the simulated genomes
		# we discard the result (these three conditions would mean that we wouldn't be studying the domain architecture, leading to bias)

		#Calculate Deletion rates ad distances within model, then push onto an array of simulated distaces
				
		my ($dels,$time,$SingleSimInterDeletionDistances) = DeletedJulianDetailed($ModelRoot,0,0,\%ModelCladeGenomesHash,$TreeCacheHash,'Simulation of DA');
		# $MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$DomArch - calculate deltion rate over tree	
					
		
		push(@$SimulatedInterDeletionDistances,@$SingleSimInterDeletionDistances);
	}
	
	#For each simulation step, calculate the observed deletion rate over the tree and the simulation steps efeectove deletion rate. Push these onto a returned distribution. Then calculate the posterior quantiles
	
	return($SimulatedInterDeletionDistances);
}

=pod
=item * RandomcalculateContinuousPosteriorQuantileModelCorrPoissonDeletionDetailed

A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.

Unlike RandomModelPoisson, however, this implementation ignores results that would produce 0 (i.e. and extinct architecture) or
100% abndance (i.e. with a deletion rate of 0 and exluded from our models).

This implementation produces more detailed output; specifically in outputs an arrayref of all the inter deletion poit distances in 

=cut

sub RandomModelCorrPoissonOptimised($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
       
    my @Weights = map{$TreeCacheHash->{$_}{'branch_length'}}@CladeGenomes; 
      
    #Construct an interval tree for fast mapping of deletion events to the tree
    #Construct an interval tree (O(nlog(n)) time) so as to computepositions in clade most efficiently 
    
    my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
    #Structure is prob => node_id
    my @Points = keys(%$ProbabilityHash);
        
    my $PointTree = Supfam::PointTree->new;
    $PointTree->build(\@Points);
  	#Create the PointTree
   
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	
	my @PoissonianDeletions = random_poisson($Iterations,$Expected_deletions); #Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my @UniformDeletions;
	
	while  (@PoissonianDeletions){ #For $Iterations number of times
	
		my $DeletionSimultation = shift(@PoissonianDeletions); #remove the first entry in array and set it as the number of deletions in this simulation
		
		@UniformDeletions = random_uniform($DeletionSimultation,0,1); # Number of deletions ($DeletionSimultation), drawn from a poissonian above, uniformly distributed across the tree.
		
		my @ModelCladeGenomes = @CladeGenomes;

		
		my $TotalDeletedGenomesHash = {};
		
		#$DeletionsNumberDistribution->{$DeletionSimultation}++; #Turn on if you want to measure how often a value is sampled
		
		foreach my $DeletionPoint (@UniformDeletions) {
                
			my $DeletionPoint = $PointTree->Search($DeletionPoint);#Find the node directly below the deletion
			my $DeletedBranch = $ProbabilityHash->{$DeletionPoint};
			#Search Tree for deletion point

			my $CurrentSimDeletedGenomes = $TreeCacheHash->{$DeletedBranch}{'Clade_Leaves'}; #Array ref to the genomes beneath the current deletion point
			
			@{$TotalDeletedGenomesHash}{@$CurrentSimDeletedGenomes} = (undef) x scalar(@$CurrentSimDeletedGenomes);
		}
		
		my $TotalDeletedGenomes = [];
		push(@$TotalDeletedGenomes,keys(%$TotalDeletedGenomesHash));
		
		## test to see if the simulation has resulted in an entire clade possessing a domain architecture and nothing else (i.e. which would make us find a new MRCA and a deletion rate of 0)
		#OR that the simulation has ended with no domain archtectures present anywhere OR that they are present everywhere, in which case we would set the deletion rate as zero
	
		my (undef,undef,undef,$ModelRemianingLeaves) = IntUnDiff($TotalDeletedGenomes,\@ModelCladeGenomes);;
	
		my $ModelFullCladeExclusive = 0; #Preallocate
		my $no_model_genomes = scalar(@$ModelRemianingLeaves);
		
		if($no_model_genomes > 0){

			my $ModelRoot = FindMRCA($TreeCacheHash,$root,$ModelRemianingLeaves);			
			(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff($TreeCacheHash->{$ModelRoot}{'Clade_Leaves'},$ModelRemianingLeaves)	; #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
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
=item * RandomModelCorrPoissonOptimised
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.

Unlike RandomModelPoisson, however, this implementation ignores results that would produce 0 (i.e. and extinct architecture) or
100% abndance (i.e. with a deletion rate of 0 and exluded from our models).

This has been optimised to run as fast as I can make it. It uses interval trees and is careful with memory. Hope that it works!
=cut

sub RandomModelCorrPoissonOptimisedDetailed($$$$$) {
	
	my ($root,$FalseNegativeRate,$Iterations,$deletion_rate,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
       
    my @Weights = map{$TreeCacheHash->{$_}{'branch_length'}}@CladeGenomes; 
      
    #Construct an interval tree for fast mapping of deletion events to the tree
    #Construct an interval tree (O(nlog(n)) time) so as to computepositions in clade most efficiently 
    
    my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
    #Structure is prob => node_id
    my @Points = keys(%$ProbabilityHash);
        
    my $PointTree = Supfam::PointTree->new;
    $PointTree->build(\@Points);
  	#Create the PointTree
   
	my $TotalBranchLength = $TreeCacheHash->{$root}{'Total_branch_lengths'};
	my $Expected_deletions = $deletion_rate*$TotalBranchLength; #$Expected_deletions is the mean of a poisson process used to model deletions
	
	my @PoissonianDeletions = random_poisson($Iterations,$Expected_deletions); #Number of deletions in this iteration. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $DeletionsNumberDistribution = {}; #This is a hash of the number of deletions modelled in the simualtion
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	my @UniformDeletions;
	
	my $SimulatedInterDeletionDistances = [];
	
	while  (@PoissonianDeletions){ #For $Iterations number of times
	
		my $DeletionSimultation = shift(@PoissonianDeletions); #remove the first entry in array and set it as the number of deletions in this simulation
		
		@UniformDeletions = random_uniform($DeletionSimultation,0,1); # Number of deletions ($DeletionSimultation), drawn from a poissonian above, uniformly distributed across the tree.
		
		my @ModelCladeGenomes = @CladeGenomes;

		
		my $TotalDeletedGenomesHash = {};
		
		#$DeletionsNumberDistribution->{$DeletionSimultation}++; #Turn on if you want to measure how often a value is sampled
		
		foreach my $DeletionPoint (@UniformDeletions) {
                
			my $DeletedBranch = $PointTree->Search($DeletionPoint);#Find the node directly below the deletion
			#Search Tree for deletion point

			my $CurrentSimDeletedGenomes = $TreeCacheHash->{$DeletedBranch}{'Clade_Leaves'}; #Array ref to the genomes beneath the current deletion point
			
			@{$TotalDeletedGenomesHash}{@$CurrentSimDeletedGenomes} = (undef) x scalar(@$CurrentSimDeletedGenomes);
		}
		
		my $TotalDeletedGenomes = [];
		push(@$TotalDeletedGenomes,keys(%$TotalDeletedGenomesHash));
		
		## test to see if the simulation has resulted in an entire clade possessing a domain architecture and nothing else (i.e. which would make us find a new MRCA and a deletion rate of 0)
		#OR that the simulation has ended with no domain archtectures present anywhere OR that they are present everywhere, in which case we would set the deletion rate as zero
	
		my (undef,undef,undef,$ModelRemianingLeaves) = IntUnDiff($TotalDeletedGenomes,\@ModelCladeGenomes);;
	
		my $ModelFullCladeExclusive = 0; #Preallocate
		my $no_model_genomes = scalar(@$ModelRemianingLeaves);
		
		my $ModelRoot;
		
		if($no_model_genomes > 0){

			$ModelRoot = FindMRCA($TreeCacheHash,$root,$ModelRemianingLeaves);			
			(undef,undef,$ModelFullCladeExclusive,undef) = IntUnDiff($TreeCacheHash->{$ModelRoot}{'Clade_Leaves'},$ModelRemianingLeaves)	; #		$ModelFullCladeExclusive will contain the members of the simulated clade beneath the simulated MRCA that aren't in the model genomes. If this is of size zero, then we should discount this result as it might incorporate bias 	
		}
		
		if ($no_model_genomes == 0  || $no_model_genomes == scalar(@CladeGenomes) || scalar(@$ModelFullCladeExclusive) == 0){

			push(@PoissonianDeletions,random_poisson(1,$Expected_deletions)); #push a number onto the end on the deletions array
			next;
		} #IFF the simulation has ended with no genomes possesing the architecture (extinction) or with complete ubiquity in the clade under study,
		# or we have ubiquity in the clade beneath the MRCA of the simulated genomes
		# we discard the result (these three conditions would mean that we wouldn't be studying the domain architecture, leading to bias)
		
		my $ModelCladeGenomesHash ={};
		map{$ModelCladeGenomesHash->{$_}=1}@$ModelRemianingLeaves;
		
		my ($dels,$time,$SingleSimInterDeletionDistances) = DeletedJulianDetailed($ModelRoot,0,0,$ModelCladeGenomesHash,$TreeCacheHash,'Simulation of DA');
		# $MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$DomArch - calculate deltion rate over tree	
					
		push(@$SimulatedInterDeletionDistances,@$SingleSimInterDeletionDistances);
		
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my ($selftest_index) =  random_uniform_integer(1,0,(scalar(@$RawResults)-1));		
	my $SelftestValue = $RawResults->[$selftest_index]; # A single uniform random simulation value
	return($SelftestValue,$distribution,$RawResults,$DeletionsNumberDistribution);
}

=pod
=item * RandomModelCorrPoissonOptimised
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.

Unlike RandomModelPoisson, however, this implementation ignores results that would produce 0 (i.e. and extinct architecture) or
100% abndance (i.e. with a deletion rate of 0 and exluded from our models).

This has been optimised to run as fast as I can make it. It uses interval trees and is careful with memory. Hope that it works!
=cut


sub calculatePosteriorQuantile($$$$){

	my ($SingleValue,$DistributionHash,$NumberOfSimulations,$CladeSize) = @_;
		
	$DistributionHash->{$SingleValue}++;
	#Stop a later step from kicking out because a there is no value in the distribution
	
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
	
	#$NumberOfSimulationsLT += ($Degeneracy/2);

	my $PosteriorQuantile = $NumberOfSimulationsLT/$NumberOfSimulations;
	
	$DistributionHash->{$SingleValue}--;
	#Undo the modification to the distribution
	
	return($PosteriorQuantile);
}


=pod * calculatePosteriorQuantile($SingleValue,$DistributionHash)

Calculates where in the total area in a probability distribution a single value occurs. Returns a value between 0 and 1. These should be uniformly distributed, from the simple fact that sum(andy distribution)
=1 and we are choosing a unifrom point from theis area.
=cut


sub calculateContinuousPosteriorQuantile($$){

	my ($SingleValue,$DistributionValues) = @_;
	
	my $NumberOfSimulationsLT = 0;
	
	map{$NumberOfSimulationsLT++ if($_ <= $SingleValue)}@$DistributionValues;
	
	my $PosteriorQuantile = $NumberOfSimulationsLT/scalar(@$DistributionValues);
	
	return($PosteriorQuantile);
}


=pod * calculateContinuousPosteriorQuantile($SingleValue,$DistributionHash)

Calculates where in the total area in a probability distribution a single value occurs. Returns a value between 0 and 1. These should be uniformly distributed, from the simple fact that sum(andy distribution)
=1 and we are choosing a unifrom point from theis area.

This is an extension from calculatePosteriorQuantile, which is a finction for discrete probability distributions
=cut

sub calculateHashContinuousPosteriorQuantile($$$){
	
	my ($SingleValue,$DistributionHash,$NumberOfSimulations) = @_;
	
	$DistributionHash->{$SingleValue}++;
	#Add the observed value ($SingleValue) to the distribution hash as a 'psuedo point'. This prevents errors kicking out because there is no point in the distribution which matches it.
	
	my $NumberOfSimulationsLT = 0;
	
	my @DistIndiciesLessThan;
	
	my $RoundedSingleVar = sprintf "%.12f", $SingleValue;

	my $flag =0;
	#A flag to ensure that a value is never seen as being equal to the single value more than once
	
	foreach my $distval (keys(%$DistributionHash)){
		
		my $rounded = sprintf "%.12f", $distval;
		
		
		print $distval." <- Not a decimal sprintf format:".$rounded."\n" unless(is_dec_number($rounded));
		#Helps with diagnosing errors
		my $dec_compare =  dec_cmp($RoundedSingleVar,$rounded);
		
		push(@DistIndiciesLessThan,$distval) if($dec_compare == -1);
		
		if($dec_compare == 0){
			
			unless($flag){
			
				$flag=1;
				
			}else{
				
				die "Two indicies in distribution are seen as being identical to the single_value. Try doing a float comparison to a greater level of precision (smaller epsilon)\n";
			}	
		}
		
		#Using dec_cmp to extract the decimal values. This prevents issues owing to floating point. Compare values to 10^-12 precision
		#To ensure that the precision level isn't too coarse grained, there's an additional step that makes sure that the single_values (observed value) isn't seen more than once.
	}
	
	
	
	foreach my $value (@DistIndiciesLessThan){
				
		$NumberOfSimulationsLT += $DistributionHash->{$value};
	}
	
	if(exists($DistributionHash->{$SingleValue})){
				
		my $Degeneracy = $DistributionHash->{$SingleValue};#Number of simulations of equal score. We place our point to sum up to uniform in this region
		my ($DegeneracyContribution) = random_uniform(1,0,$Degeneracy);
				
		$NumberOfSimulationsLT += $DegeneracyContribution;
		#$NumberOfSimulationsLT += ($Degeneracy/2);
	}
		
	my $PosteriorQuantile = $NumberOfSimulationsLT/($NumberOfSimulations+1);
	#Add plus one to the number of simulations as we added a further 

	die "Posterior Quantile > 1!! This shoudl never ever happen!  Post Quant: $PosteriorQuantile, Nsims = $NumberOfSimulations, NsimsLT = $NumberOfSimulationsLT" if ($PosteriorQuantile > 1);
	
	$DistributionHash->{$SingleValue}--;
	#Remove the previosly added psuedo value.
	
	return($PosteriorQuantile);
}

=pod * calculateContinuousPosteriorQuantile($SingleValue,$DistributionHash)

Calculates where in the total area in a probability distribution a single value occurs. Returns a value between 0 and 1. These should be uniformly distributed, from the simple fact that sum(andy distribution)
=1 and we are choosing a unifrom point from this area.

This is an extension from calculatePosteriorQuantile, which is a finction for discrete probability distributions

This function uses a hash to loop through the distribution. Pass this into the funtion
=cut



sub calculateOldStylePosteriorQuantile($$$$){

	my ($SingleValue,$distribution,$NumberOfSimulations,$CladeSize) = @_;
		
  my $CumulativeCount= 0;
  my @CumulativeDistribution; #(P(Nm<nr)) nr (number of genomes in reality) is the index, Nm is the random variable

	
#Create cumulative distibution list
      for my $NumberOfGenomes (0 .. $CladeSize){
		$CumulativeCount += $distribution->{$NumberOfGenomes} if(exists($distribution->{$NumberOfGenomes})); # Cumulative count is a sum of the frequency of genome observation up to this point	
        $CumulativeDistribution[$NumberOfGenomes]=$CumulativeCount;  
 	
      }
 	
  	
    my $ProbLT = ($CumulativeDistribution[$SingleValue-1])/$NumberOfSimulations;#Cumlative probability of less genomes in model simulations than in reality
    
	return($ProbLT);
}


=pod * calculateOldStylePosteriorQuantile($SingleValue,$DistributionHash)

An older implementation of the posterior quantile test that resulted in incorrect values. Kept in the repo for legacy.

=cut

sub calculateJulianStylePosteriorQuantile($$$$){

	my ($SingleValue,$dist,$NumberOfSimulations,$CladeSize) = @_;
		
	
my %distribution=%$dist;

my $ii=0;
my $gennum = $SingleValue;


for my $i (0 .. $CladeSize){
	
if ($i == $gennum){

	last;
}

if (exists($distribution{$i})){

	$ii=$ii+$distribution{$i};
  }
}

	return($ii/$NumberOfSimulations);
}


=pod * calculateJulianStylePosteriorQuantile($SingleValue,$DistributionHash)

An older implementation of the posterior quantile test as per Julian Gough's orginal implementation of hgt.pl. Should be equivilent to calculateOldStylePosteriorQuantile

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
