package Supfam::hgtoptimised;
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
				Clade 
				Deleted
				RandomModel
				RandomModelJulian
				DeletedJulian
				DeletedPoisson
				RandomModelPoisson
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;
use Bio::TreeIO;

use Bio::Tree::TreeFunctionsI;
use Time::HiRes;
use Math::Random;
use Math::Random qw(random_poisson);

sub Clade($$) {
	
	my ($tree,$CladeGenomes) = @_;
	
	my $flag = 1;
	#This routine includes an iterative process to find the common ancestor of all genomes in @$CladeGenomes. This is the flag
	
	my $AncestorNodeID = $CladeGenomes->[0]; # An initial value 
	
	
		while($flag){
			
		my $node = $tree->find_node(-id => "$AncestorNodeID");
		my @Descendents = $node->get_all_Descendents;
		
		IntUnDiff($CladeGenomes,\@Descendents);		
		my $Outgroup = $_[2];
		#This is a pointer to an array of elements unique to the first input list to IntUnDiff (so @$CladeGenomes).
		
			if(scalar(@$Outgroup)){
				
			$AncestorNodeID = $node->ancestor;
				
			}else{
				$flag =0;
			}
		}

return($AncestorNodeID);

}


=pod
=item * Clade($$)
Given a tree and a list of genomes, this function will find the most recent common ancestor of those genomes in the tree given. Clade($tree,$PointerToGenomesArray)
=cut

sub Deleted{ #Left out the prototyping arguments as its a recursive function. Be careful
	
	my ($tree,$node,$dels,$time,$GenomesOfDomArch) = @_;
	#$time in this function is evolutionary time and $dels are the number of observed deletions in that time
	
	my $root = $tree->get_root_node;
	
	my @Descendents = map{$_->id} grep{ $_->is_Leaf==1 } $node->get_all_Descendents ;# All descendent leave ids. Sorry about it looking ugly, but this was the easiest way, honest!
	@Descendents = ($node->id) if(($node->is_Leaf)); # This allows for when the node is a leaf and hence has no descendents
	
	my $UnionHash={};
	my $Intersection =[];

	foreach my $element (@Descendents, @$GenomesOfDomArch) { 
		if($UnionHash->{$element}++ == 2){
			push(@$Intersection,$element);
		} 
	}# Note that this method requires that the two lists of genomes be unique (no internal repeats in the lists)!

	if(scalar(@$Intersection)){ # if there is any overlap between the list of genomes containing dom arch and those in the clade being studied
	
		my $EvolutionaryTime = 0; # Evolutionary time from ancestor
		$EvolutionaryTime = $node->branch_length unless($node eq $root);
		
		$time += $EvolutionaryTime;
	
		unless ($node->is_Leaf){
			
			my @Children = $node->each_Descendent; # All the immediate children of the node
			
				foreach my $Child (@Children){
			
					($dels,$time) = Deleted($tree,$Child,$dels,$time,$GenomesOfDomArch);
					#Notice the replacement of $node with $Child.
				}
				
			return($dels,$time);
				
		}else{
			
			#We haven't seen a deletion. Don't increase $dels and return
			return($dels,$time);	
		}
		
	}else{
		
		my $EvolutionaryTime = $node->branch_length; # Evolutionary time from ancestor
		
		$time += $EvolutionaryTime/2; #Assume that the deletion took place at some point between last ancestor and now.
		$dels++;
		#We've seen a deletion. Increase $dels, update $time and return
		
		return($dels,$time);
	}
}


=pod
=item * Deleted
This function looks at deletions of a particular domain archtecture of interest. Given a tree and a list of genomes that posses that domain architecture, this function
moves through the tree, level by level, clade by clade. If a clade contains an example of the domain arch, add the evolutionary time distance of the clade root from its parent
and look at the clades of its children nodes. If no examples of the dom arch are found, increase the number of deletions by one and move on. 
Deleted($BioPerltreeobject,$CladeParentNode,$dels,$time,$GenomesPossessingDomainArchitecture)
=cut

sub DeletedJulian($$$$$$$){ 
	#($subtree,0,0,$NodesObserved)
	my ($MRCA,$dels,$time,$HashOfGenomesObserved,$TreeCacheHash,$treeroot,$DomArch) = @_;
	#$time in this function is evolutionary time and $dels are the number of observed deletions in that time
		
    my @InternalNodes = @{$TreeCacheHash->{$MRCA}{'each_Descendent'}};
    @InternalNodes = $MRCA if($TreeCacheHash->{$MRCA}{'is_Leaf'});
    #InternalNodes will be a list of all nodes within the tree, grown in the loop below. Some initial values are thrown in here.
           
	    while(scalar(@InternalNodes)){
	    	
	    	my $node = pop(@InternalNodes);
	    	
	    	my $DistanceFromAncestor = $TreeCacheHash->{$node}{'branch_length'};
 	
	    	my @Descendents = map{$_->id} @{$TreeCacheHash->{$node}{'Clade_Leaves'}}  ;# All descendent leave ids. Sorry about it looking ugly, but this was the easiest way, honest!
			@Descendents = ($node->id) if($TreeCacheHash->{$node}{'is_Leaf'}); # This allows for when the node is a leaf and hence has no descendents
	  	   
			my $Intersection =[];
			no warnings 'uninitialized';
			foreach my $descendent (@Descendents) { 
				if($HashOfGenomesObserved->{$descendent} == 1){
					push(@$Intersection,$descendent);
				}elsif($HashOfGenomesObserved->{$descendent} > 1){
					print $DomArch."\n";
					die '$HashOfGenomesObserved should have UNIQUE key value pairs. '.$descendent.' is non unique. There is a serious issue with the way that data is inputed into the hash';
				}
			}# Note that this method requires that the two lists of genomes be unique (no internal repeats in the lists)!
	   	   	
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
	    
	    #Consider returning upper and lower estimation confidence intervals
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
	
	my ($root,$FalseNegativeRate,$Iterations,$NoDeletions,$TreeCacheHash) = @_;
	
	#$root is the root of the subtree or the most recent common ancestor

    my @CladeGenomes = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
    push(@CladeGenomes,$root) if ($TreeCacheHash->{$root}{'is_Leaf'});
    my %CladeGenomesHash;
    map{$CladeGenomesHash{$_} =1;}@CladeGenomes ;#Initialise a hash of the genomes in this subtree
           
	my $ProbabilityHash = $TreeCacheHash->{$root}{'Probability_Hash'};
	my @ProbabilityIntervals = sort(keys(%$ProbabilityHash));
	
	my @PoissonianDeletions = random_poisson($Iterations,$NoDeletions); #Number of deletions. This is drawn from a poissonian with mean equal to the number of deletions (the MLE)
	
	my $RawResults = []; #Create an array to store the direct simulation results, rather than the results aggregated into a hash like $distribution  
	my $distribution = {}; # This is ultimately what the distributon of the model runs will be stored in
	
	for my $DeletionSimultation (@PoissonianDeletions){ #For $Iterations
		
		my @UniformDeletions = random_uniform($DeletionSimultation); # Number of deletions, drawn from a poissonian above, uniformly distributed across the tree
		my %ModelCladeGenomesHash = %CladeGenomesHash;
		
		foreach my $DeletionPoint (@UniformDeletions) {
                    
			my $index = 0; while ($DeletionPoint > $ProbabilityIntervals[$index]){$index++;}
			# @ProbabilityIntervals is a precalculated hash of all the nodes in the sub-tree from the MRCA ($root) and where they sit in a stretched out sum of all branch lengths. $DeletedNode is a uniform random point along this line
			my $DeletedNode = $ProbabilityHash->{$ProbabilityIntervals[$index]};
			
			map{delete($ModelCladeGenomesHash{$_})}@{$TreeCacheHash->{$DeletedNode}{'Clade_Leaves'}};
			delete($ModelCladeGenomesHash{$DeletedNode}) if ($TreeCacheHash->{$DeletedNode}{'is_Leaf'});		
		}
		
		my $no_model_genomes = scalar(keys(%ModelCladeGenomesHash));
		
		$distribution->{$no_model_genomes}++;
		push(@$RawResults,$no_model_genomes);
		#Update the distribution of the run accordingly and store results in rawresults
	}
	
	my $SelftestValue = $RawResults[scalar(rand(@$RawResults))]; # A single uniform random simulation value
		
	return($SelftestValue,$distribution,$RawResults);
}

=pod
=item * RandomModelPoisson
A deletion model based on the poisson distribution of deletion events. Using the number of deletions in the entire clade to parameterise the distribution (as obtainined using DeletedJulian),
we draw $itr intergers from the distribution and then scatter then, for each of those numbers, scatter N deletion events over the tree, where N is the poisson number.
=cut

sub RandomModel {
	
	my ($CladeGenomes,$ObservedGenomes,$TreeCacheHash,$DeletionRate,$MRCA) = @_;
	
	#Note, these $*Genome variables are BioPerl nodeID objects
	
	my $LHHash={};
	$LHHash->{$MRCA}=1;
	
	my $TreeLikelihood = 1;
	
	#print STDERR "Observed genome likelihoods:\n";
	
#	my $Observed = join("\n",@$ObservedGenomes);
#	open OBS, ">ObservedID.dat";
#	print OBS $Observed;
#	close OBS;
#	

	my @UnobservedIDs;

	foreach my $CladeGenomeID (@$CladeGenomes){
		
		if(grep{/$CladeGenomeID/}@$ObservedGenomes){ #If $CladeGenome is in the list of genomes containing architecture
			
		my $ObservedGenomeLH = PositiveDomArchObservations($LHHash,$CladeGenomeID,$TreeCacheHash,$DeletionRate); #Populate $LHHash with LH values resulting from observed values
		$TreeLikelihood = $TreeLikelihood*$ObservedGenomeLH;
		
		#my $NodeName = $CladeGenomeID->id; # Retrieve the real (SUPERFAMILY) name of the leaf in question	
		#print STDERR $NodeName.' : '.$ObservedGenomeLH."\n";
			
		}else{
			
			push(@UnobservedIDs,$CladeGenomeID);
			}
	}

	#print STDERR "Unobserved genome likelihoods:\n";
	
	foreach my $UnobservedGenome  (@UnobservedIDs){

		my $UnobservedGenomeLH = NegativeDomArchObservations($LHHash,$UnobservedGenome,$TreeCacheHash,$DeletionRate); #Populate $LHHash with LH values
		$TreeLikelihood = $TreeLikelihood*(1-$UnobservedGenomeLH); #The value is 1-prob as we are dealing with the probability of NOT observing dom arch
		
		#my $NodeName = $UnobservedGenome->id;
		#print STDERR $NodeName.' : '.$UnobservedGenomeLH."\n";
	}
	
	return($TreeLikelihood);

}

=pod
=head2 Methods
=over 4
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
be one, as it is the only possible way to explain the previously observed sighiting of the do arch.
=cut

1;