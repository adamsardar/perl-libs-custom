package Supfam::TreeFuncs;
require Exporter;

=head1 NAME

Supfam::TreeFuncs

=head1 SYNOPSIS

Holds functions related to parsing Newick trees and calculating domain architecture information 
use Supfam::TreeFuncs;

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
			CalculateLineage
			ExtractSubtree
			FindTrueRoot
			BuildTreeCacheHash
			GenerateCladeTimeHash
			SQL2Newick
			TreeIntersection
			TreeHash2Newick
			FindMRCA
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use DBI;
#use Data::Dumper;
#use Term::ProgressBar;
#use Math::Combinatorics;

use Supfam::Config;
use Supfam::Utils;
use Supfam::SQLFunc;

use Bio::TreeIO;
use Bio::Tree::TreeFunctionsI;
use IO::String;
use Time::HiRes;
use Bio::Tree::TreeI;


=pod
=head2 Methods
=over 4
=cut


sub SQL2Newick($$){
		
	my ($rootleft,$rootright) = @_;
		
	my $RightTreeHash = {}; #A hash of the tree ordered by right id
	my $LeftTreeHash = {}; #A hash of the tree ordered by left id
	
	my $dbh = dbConnect();
	
	my $sth = $dbh->prepare("SELECT left_id,right_id,edge_length,nodename FROM tree WHERE left_id >= ? AND right_id <= ?;");
	$sth->execute($rootleft,$rootright);
	
	while (my ($leftid,$rightid,$edge_length,$nodename) = $sth->fetchrow_array() ){
	
		my $NewickClade = 'NULL'; #The calde below this node in newick format - initialise with a default value
		
		my $NodeData = [$nodename,$edge_length,$NewickClade];
		
		$RightTreeHash->{$rightid} = [$leftid,$NodeData];
		$LeftTreeHash->{$leftid} = [$rightid,$NodeData];
	}
	
	my $NextGenLeftIDs = {};
	my $CurrentGenLeftIDs = [];
	
	my $TreeLeftIDs; 
	map{$TreeLeftIDs->{$_}=undef}keys(%$LeftTreeHash);
	
	foreach my $leftid (keys(%$TreeLeftIDs)){
		
		my $rightid = $LeftTreeHash->{$leftid}[0];
		
		if ($rightid == ($leftid+1)){ #if node is a leaf
				
			no warnings 'uninitialized';	
				
			my ($nodename,$edge_length,$NewickClade) = @{$LeftTreeHash->{$leftid}[1]};
				
			my $NewickString = $nodename.":".$edge_length;
			$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string from 'NULL' to correct value
			$RightTreeHash->{$rightid}[1][2] = $NewickString;
		
			delete($TreeLeftIDs->{$leftid});
		}	
	}
	
	#Initilise leaf nodes with values
	#(possibility of additional functionality - exclude a list of genomes)
		
	while (scalar(keys(%$TreeLeftIDs))){

	
		
	my ($rootleft,$rootright) = @_;
		
	my $RightTreeHash = {}; #A hash of the tree ordered by right id
	my $LeftTreeHash = {}; #A hash of the tree ordered by left id
	
	my $dbh = dbConnect();
	
	my $sth = $dbh->prepare("SELECT left_id,right_id,edge_length,nodename FROM tree WHERE left_id >= ? AND right_id <= ?;");
	$sth->execute($rootleft,$rootright);
	
	while (my ($leftid,$rightid,$edge_length,$nodename) = $sth->fetchrow_array() ){
	
		my $NewickClade = 'NULL'; #The calde below this node in newick format - initialise with a default value
		
		my $NodeData = [$nodename,$edge_length,$NewickClade];
		
		$RightTreeHash->{$rightid} = [$leftid,$NodeData];
		$LeftTreeHash->{$leftid} = [$rightid,$NodeData];
	}
	
	my $NextGenLeftIDs = {};
	my $CurrentGenLeftIDs = [];
	
	my $TreeLeftIDs; 
	map{$TreeLeftIDs->{$_}=undef}keys(%$LeftTreeHash);
	
	foreach my $leftid (keys(%$TreeLeftIDs)){
		
		my $rightid = $LeftTreeHash->{$leftid}[0];
		
		if ($rightid == ($leftid+1)){ #if node is a leaf
				
			no warnings 'uninitialized';	
				
			my ($nodename,$edge_length,$NewickClade) = @{$LeftTreeHash->{$leftid}[1]};
				
			my $NewickString = $nodename.":".$edge_length;
			$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string from 'NULL' to correct value
			$RightTreeHash->{$rightid}[1][2] = $NewickString;
		
			delete($TreeLeftIDs->{$leftid});
		}	
	}
	
	#Initilise leaf nodes with values
	#(possibility of additional functionality - exclude a list of genomes)
		
	while (scalar(keys(%$TreeLeftIDs))){

		foreach my $leftid (keys(%$TreeLeftIDs)) {
			
			my $rightid = $LeftTreeHash->{$leftid}[0];
			
			my ($LeftNewick,$RightNewick) = ('NULL','NULL');
			$LeftNewick = $LeftTreeHash->{$leftid+1}[1][2] if (exists($LeftTreeHash->{$leftid+1}));
			$RightNewick = $RightTreeHash->{$rightid-1}[1][2] if (exists($RightTreeHash->{$rightid-1}));
			#Extract the newick format trees of the nodes to the left and right descending edges away from this node
									
			if($LeftNewick ne 'NULL' && $RightNewick ne 'NULL'){
				
				#Aggregate two nodes below this one into a newick format
				my $edge_length = $LeftTreeHash->{$leftid}[1][1];
				my $NewickString = "(".$LeftNewick.",".$RightNewick."):".$edge_length;
				$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string
				$RightTreeHash->{$rightid}[1][2] = $NewickString;
							
			}
					
			delete($TreeLeftIDs->{$leftid}) if($LeftTreeHash->{$leftid}[1][2] ne 'NULL');	
		}	
	}
	
	##Ugly!
	my $FullTree = $LeftTreeHash->{$rootleft}[1][2];
	$FullTree="(".$FullTree.");";
	$LeftTreeHash->{$rootleft}[1][2]=$FullTree;
	$FullTree = $RightTreeHash->{$rootright}[1][2];
	$FullTree="(".$FullTree.");";
	$RightTreeHash->{$rootright}[1][2]=$FullTree;
	#Above is to correct the final trees created in the resulting hashes.
	
	return($RightTreeHash,$LeftTreeHash);
		foreach my $leftid (keys(%$TreeLeftIDs)) {
			
			my $rightid = $LeftTreeHash->{$leftid}[0];
			
			my ($LeftNewick,$RightNewick) = ('NULL','NULL');
			$LeftNewick = $LeftTreeHash->{$leftid+1}[1][2] if (exists($LeftTreeHash->{$leftid+1}));
			$RightNewick = $RightTreeHash->{$rightid-1}[1][2] if (exists($RightTreeHash->{$rightid-1}));
			#Extract the newick format trees of the nodes to the left and right descending edges away from this node
									
			if($LeftNewick ne 'NULL' && $RightNewick ne 'NULL'){
				
				#Aggregate two nodes below this one into a newick format
				my $edge_length = $LeftTreeHash->{$leftid}[1][1];
				my $NewickString = "(".$LeftNewick.",".$RightNewick."):".$edge_length;
				$LeftTreeHash->{$leftid}[1][2] = $NewickString; #Update old newick string
				$RightTreeHash->{$rightid}[1][2] = $NewickString;
							
			}
					
			delete($TreeLeftIDs->{$leftid}) if($LeftTreeHash->{$leftid}[1][2] ne 'NULL');	
		}	
	}
	
	##Ugly!
	my $FullTree = $LeftTreeHash->{$rootleft}[1][2];
	$FullTree="(".$FullTree.");";
	$LeftTreeHash->{$rootleft}[1][2]=$FullTree;
	$FullTree = $RightTreeHash->{$rootright}[1][2];
	$FullTree="(".$FullTree.");";
	$RightTreeHash->{$rootright}[1][2]=$FullTree;
	#Above is to correct the final trees created in the resulting hashes.
	
	return($RightTreeHash,$LeftTreeHash);
}

=pod
=item *SQL2Newick(root left id, root right id)

Given the left_id and the right_id of a node, this will return the tree stored in SUPERFAMILY under that node in the tree table as a newick string.
Both right and left ids need to be specified so as to prevent the all to easy error choosing the wrong node as the root. Returns two hashes: one has
keys of left_ids and which hash to: $LeftHash->{leftid}=[rightid,[($SQLnodename,$edge_length,$CladeAsNewickString)]]. The other hash is the same but with 
rightids.

This function only deals with binary trees, as in SUPERFAMIY.
=cut

sub BuildTreeCacheHash($){
	
	my ($TreeFile) = @_;
	
	my $tic = Time::HiRes::time; 
	
	my $TreeCacheHash = {}; # A massive limitation of this script is the use of BioPerl TreeIO which is SLOOOOW. This is a lookup hash to speed things up.
	
	my $input = new Bio::TreeIO(-file   => "$TreeFile",
	                            -format => "newick") or die $!;
	                            
	my $tree = $input->next_tree;
	my $root = $tree->get_root_node;
	
	#Read in and initialise tree
	
	my $BranchLengthsFlag = 0; #Flag to test if the input file contained branch length info (will no compute cladetime hash if not)	
	
	## Initialise all the tree values that we might need in this hash to minimise calls to BioTree.
	
	$TreeCacheHash->{$root}={};
	my @RootDescendents = $root->get_all_Descendents;
	$TreeCacheHash->{$root}{'all_Descendents'}=\@RootDescendents;
	map{$TreeCacheHash->{$_}={}}@{$TreeCacheHash->{$root}{'all_Descendents'}};
	#Bootstrap the BioPerl internal nodeids as keys to Cache Hash
	
	foreach my $node ($root->get_all_Descendents,$root) {
	
		#each_descendant
		$TreeCacheHash->{$node}{'each_Descendent'} = [];
		push(@{$TreeCacheHash->{$node}{'each_Descendent'}},$node->each_Descendent);
		
		#alldescendents
		my @AllDescendents = $node->get_all_Descendents;
		$TreeCacheHash->{$node}{'all_Descendents'} =  \@AllDescendents;
		
		#ancestor
		my $Ancestor = $node->ancestor;
		$TreeCacheHash->{$node}{'ancestor'} =  $Ancestor unless($node eq $root);
		
		#nodeid (as in the tree file provided - Real World ID)
		my $RWnodeID = $node->id;
		$TreeCacheHash->{$node}{'node_id'} =  $RWnodeID; 
		
		#Clade_leaves
		my @CladeLeaves = grep{$_->is_Leaf == 1}@{$TreeCacheHash->{$node}{'all_Descendents'}};
		$TreeCacheHash->{$node}{'Clade_Leaves'} = \@CladeLeaves;
		
		#Branch Length from ancestor from ancestor
		my $Branch = $node->branch_length;
		
		$TreeCacheHash->{$node}{'branch_length'} = 0 if ($node eq $root);
		
		print STDERR "Node id:".$TreeCacheHash->{$node}{'node_id'}." - branch length of zero found, introducing a pseudocount of 0.000001\n"  if ($Branch ~~ 0 && $node ne $root);				
		$Branch = 0.000001 if ($Branch ~~ 0 && $node ne $root); 
		
		$TreeCacheHash->{$node}{'branch_length'} = $Branch ;
		
		$BranchLengthsFlag = 1 unless($Branch ~~ 0 || $Branch ~~ undef);
		
		#Total Branch Lengths in tree - I'd love a section in here: maybe in the future
	}
	
	#is_Leaf Entries
	map{$TreeCacheHash->{$_}{'is_Leaf'} = 0}@{$TreeCacheHash->{$root}{'all_Descendents'}};
	map{$TreeCacheHash->{$_}{'is_Leaf'} = 1}@{$TreeCacheHash->{$root}{'Clade_Leaves'}};

	if($BranchLengthsFlag){
	
		map{GenerateCladeTimeHash($_,$TreeCacheHash)}(@{$TreeCacheHash->{$root}{'all_Descendents'}},$root);
	}
	#This adds two other key/value pairs:
	# 	$TreeCacheHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
	# 	$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
	
	#$BranchLengthsFlag is only for the case where the input tree has no branch lengths;
	
	my $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Tree Cache hash:".($toc-$tic)."seconds\n";
	
	return($root,$TreeCacheHash, $tree);
}

=pod
=item * BuildTreeCacheHash($TreeFile)
Often it is much simpler to simply provide a string of a newick tree to this subroutine and recieve a hash of that tree with a whole load
of information about each node precalculated. Speed is one possible reason to do this in case there are many calls to bioperl occuring in your script.
There are also useful values that are not diectly given via bioperl.

Each node is stored as follows $TreeHash->{BioPerl NodeID}
There is then a variety of values stored:
$TreeHash->{BioPerl NodeID}{'all_Descendents'} - All BioPerl Node IDs of nodes below this one
$TreeHash->{BioPerl NodeID}{'each_Descendent'} - The direct descendents of this node
$TreeHash->{BioPerl NodeID}{'Clade_Leaves'} - All the leaf nodes below this one
$TreeHash->{BioPerl NodeID}{'branch_length'} - The length of the branch connecting this node to it's parent
$TreeHash->{BioPerl NodeID}{'is_Leaf'} - 1/0 flag for if the node is a leaf
$TreeHash->{BioPerl NodeID}{'ancestor'} - BioPerl reference of nodes parent
$TreeCacheHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
$TreeCacheHash->{$MRCA}{'node_id'}='Node name' as given in the input newick file

returned is:

($root,$TreeCacheHash, $tree);

$root - the bioperl ID of the root

$TreeCacheHash - as above

$tree - the BioPerl representation of the tree

=cut

sub GenerateCladeTimeHash($$){
	
	#Provide a node and the tree hash, return a pointer to the hash of {cladepoint => node}. Cladepoint is equal to the normalised (0<x<1) point in the total clade time at which the node sits
	my ($MRCA, $TreeCacheHash) = @_;
	#MRCA is the most recent common acnestor of clade under study

	#Total Branch Length in clade below this point
	my $TotalBranch = 0;
	map{$TotalBranch += $TreeCacheHash->{$_}{'branch_length'}}@{$TreeCacheHash->{$MRCA}{'all_Descendents'}}; #Sum all branch lengths beneath MRCA of clade
	$TreeCacheHash->{$MRCA}{'Total_branch_lengths'} = $TotalBranch ;
	$TreeCacheHash->{$MRCA}{'Total_branch_lengths'} = 0 if ($TreeCacheHash->{$MRCA}{'is_Leaf'});
		
	my $time = 0;
	my $ProbabilityIntervalHash = {};

	unless($TreeCacheHash->{$MRCA}{'is_Leaf'}){
			
		my @nodes = @{$TreeCacheHash->{$MRCA}{'each_Descendent'}};
		
		while (my $node = pop(@nodes)){
					
				@nodes=(@nodes,@{$TreeCacheHash->{$node}{'each_Descendent'}}); 
				# Add the children of the internal node under study in this loop to the list, unless it's a leaf node
	
				$time += $TreeCacheHash->{$node}{'branch_length'};
				
				#Sometimes a tree can have a branch length of zero - particularly when the tree has been renormalised for deletions

				$ProbabilityIntervalHash->{$time/$TotalBranch} = $node;
		}
	}else{
		
		$ProbabilityIntervalHash = 'NULL';
	}
	
	$TreeCacheHash->{$MRCA}{'Probability_Hash'} = $ProbabilityIntervalHash;
		
	return(1);
}

=pod
=item * GenerateCladeTimeHash($MRCA, $TreeCacheHash)
An internal function for BuildTreeCacheHash. Adds two additional values onto hash:
$TreeCacheHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
Don't call this function directly please.No packages will be installed, upgraded, or removed.
0 packages upgraded, 0 newly installed, 0 to remove and 52 not upgraded.

=cut

sub CalculateLineage($$$$); #Predeclare so as to stop perl complainign about being 'called too early to check prototype'

sub CalculateLineage($$$$) {
	
	my ($Lineage,$TaxID,$LCA, $dbh) = @_;
	$dbh = dbConnect() unless defined $dbh;
	my $close_dbh = (@_ < 3)?1:0;
	
	my $query = $dbh->prepare_cached("select parent_id,name,rank  from ncbi_taxonomy where taxon_id = ?") or die "Query failed";
	$query->execute($TaxID) or return undef;
	
	my @SQLOutput = $query->fetchrow_array();
	$query -> finish;
	
	my ($parent,$name,$rank) = @SQLOutput;
	$Lineage = $name.'::'.$Lineage;
	
		if($TaxID == $LCA){
			 
			return ($Lineage);
	
		}elsif($TaxID == 0){
			
			return $Lineage;
			print STDERR "Warning, Last Common Ancestor provided was not found before reaching the tree route. The code might not be working as expected";
		}else{
		
			my $output = CalculateLineage($Lineage,$parent,$LCA, $dbh);
			return ($output);
		}
		
	dbDisconnect($dbh) if $close_dbh;
}

=pod
=item * CalculateLineage($Lineage,$TaxID,$LCA, $dbh)
A recursive funtion that aims to find the full NCBI lineage of a NCBI taxon, given by $TaxID. $dbh is optional (a connection to the SUPERFAMILY database using dbconnect will be made). When called,
$Lineage can be set as blank ('') as it is populated through recursive callign of the script
=cut

sub ExtractSubtree($$$){
	
	my ($tree,$Ingroup,$remove) = @_; 
	#Parse in a bioperl tree object reference, a refernce to an array of the ingroup taxa that you wish to have in resulting tree, a flag as to whehter to
	#include other taxa beneath the ingroup MRCA (0) or not (1).

	my $root = $tree->get_root_node;
	
	my @TreeLeafNodeIDs = grep{$_->is_Leaf}$root->get_all_Descendents;
	my $IngroupNodeIDs = [];
	map{my $NODENAME = $_->id; push(@$IngroupNodeIDs,$_) if grep{m/^$NODENAME$/}@$Ingroup;}@TreeLeafNodeIDs; #Extract the BP node-ids for the ingroup
	
	print scalar(@$IngroupNodeIDs)."\n";
	my $MRCA = $tree->get_lca(-nodes => $IngroupNodeIDs) ; #Find most Recent Common Ancestor
	
	my $subtree = Bio::Tree::Tree->new(-root => $MRCA, -nodelete => 1);
	
	my $NewRoot = $subtree->get_root_node;
	my @NewTreeLeafNodeIDs = grep{$_->is_Leaf}$NewRoot->get_all_Descendents;
	my @NewTreeLeafNames = map{$_->id}@NewTreeLeafNodeIDs;

	$IngroupNodeIDs = [];
	map{my $NODENAME = $_->id; push(@$IngroupNodeIDs,$_) if grep{m/^$NODENAME$/}@$Ingroup;}@NewTreeLeafNodeIDs;
	my @PresentIngroup = map{$_->id}@$IngroupNodeIDs;

	if($remove){
	
		$subtree->splice(-remove_id => \@NewTreeLeafNames, -preserve_lengths => 1,-keep_id => \@PresentIngroup);
		#Remove all nodes that weren't directly specified by the ingroup
	
		#Tidy Up Dangling leaves (i.e. those with no name)
		#Delete leaf nodes with no name (i.e. internal nodes that no-longer have any children and are hence pointless) and follow them up the tree, remove further leaf nodes with no id
		
		my @DanglingNodes = grep{$_->is_Leaf}$subtree->find_node(-id => '');

		while (scalar(@DanglingNodes) > 0 ){
	
			foreach my $DanglingNode (@DanglingNodes){
		
				$subtree->remove_Node($DanglingNode);
			}
			
			@DanglingNodes = grep{$_->is_Leaf}$subtree->find_node(-id => '');
		}
	
		#Trim out useless nodes with only one child (i.e connect grandparent to it's grandchild if the intervening partent node has only one ancestor)
		my @PointlessNodes = grep{! $_->is_Leaf}grep{scalar($_->each_Descendent) == 1}($NewRoot->get_all_Descendents);
		map{$_->id('Pointless')}@PointlessNodes;
		$subtree->splice(-remove_id => 'Pointless', -preserve_lengths => 1); #Frustratingly, BioPerl seems to want to work by node names here, which are non-unique. Hmpf.
	}
	
return($subtree);
}

=pod
=item * ExtractSubtree(BioPerl TreeObject, Ingroup Array Pointer, Flag As To Delete All Non-Ingroup Genomes From Ouput Tree)
Given a tree and an ingroup, extract a subtree with root as the MRCA of the ingroup.Returns a bioperl tree object of the subtree.
If the third argument ($remove) is set to 1, then only genomes in the ingroup will appear in the returned tree, else all genomes in the original tree below the MRCA will be returned.
=cut

sub FindTrueRoot($$) {
	
	my ($tree,$OutgroupNodeIDs) = @_;
	
	$tree->reroot($$OutgroupNodeIDs[rand(scalar(@$OutgroupNodeIDs))]); #Reroot trree on a random outgroup node. It doesn't really matter which one
	
	my $TreeLeaves = [($tree->get_leaf_nodes)];
	my ($Union,$Intersection,$Ingroup,$OutgroupExclusive) = IntUnDiff($TreeLeaves,$OutgroupNodeIDs);
	
#	die "Outgroup not wholly contained within tree taxa!\n" if (scalar(@$OutgroupExclusive)); # $OutgroupExclusive shoudl be empty!
	
	my $ExitFlag = 0; #While loop below will exit when a node is found with only ingroup descndents (note that descendents has a wierd meaning here as we have rooted the tree on a leaf)
	
	my $root = $tree->get_root_node;
	my @rootdescendents = $root->each_Descendent;
	my $CurrentGenParent = $rootdescendents[0]; #As this is rooted on a leaf, there is only one descendent
	my ($NextGenParent,@CurrentGenChildren);
	
	#Test if input tree is rooted or unrooted: if rooted, then node will have exactly two descendents (as it's rooted at a midpoint). If unrooted, then it will have more or 1 (if it's a leaf)
	
	my $RootedTest = (scalar(@rootdescendents) == 2)?1:0;
	
	while (! $ExitFlag){

		@CurrentGenChildren = ($CurrentGenParent->each_Descendent);
		
		### Calculate the overlap of two pairs of two sets - descndants of the first child of this parent node and the ingroup and the decendents of the second child ...		
						
		my $FirstChildLeafDescendents = [grep{$_->is_Leaf}(($CurrentGenChildren[0])->get_all_Descendents)];
		my ($UnionFirstIngroup,$IntersectionFirstIngroup,$IngroupExclusiveFirst,$FirstChildExclusive) = IntUnDiff($Ingroup,$FirstChildLeafDescendents);
		
		my $SecondChildLeafDescendents = [grep{$_->is_Leaf}(($CurrentGenChildren[1])->get_all_Descendents)];
		my ($UnionSecondIngroup,$IntersectionSecondIngroup,$IngroupExclusiveSecond,$SecondChildExclusive) = IntUnDiff($Ingroup,$SecondChildLeafDescendents);
		
		###
		my $NextGenDesc;
		## Choose which child node to use as next gen parent (we are moving down theough the tree using a geedy algorithm)
		
		if(scalar(@$IntersectionFirstIngroup) > scalar(@$IntersectionSecondIngroup)){
			
			$NextGenParent = $CurrentGenChildren[0];
			$NextGenDesc = $FirstChildLeafDescendents;
			
		}elsif(scalar(@$IntersectionFirstIngroup) < scalar(@$IntersectionSecondIngroup)){
			
			$NextGenParent = $CurrentGenChildren[1];
			$NextGenDesc = $SecondChildLeafDescendents;
		}else{
			
			die "Poor choice of outgroup! At one node ($CurrentGenParent) there are qually as many ingroup nodes in either child\n";
		}	
		#Test for overlap of next gen descendents with the outgroup - if there are no members then we have found our root.
		
		$CurrentGenParent = $NextGenParent;
		
		my ($UnionNextGenOutgroup,$IntersectionNextGenOutgroup,$NextGenExclusive,$OutgroupExclusive) = IntUnDiff($NextGenDesc,$OutgroupNodeIDs);
		$ExitFlag = 1 unless(scalar(@$IntersectionNextGenOutgroup)); #If there are no members of the out group in the descendents list, we have our root! 
		
	}
	
	my $OldRoot = $root;
	
	print map{$_->id}($NextGenParent->each_Descendent);
	print "\n";
	print $root->id;
	print "\n";
	
	$tree->reroot_at_midpoint($NextGenParent);
	
	$tree->splice(-remove_id => [$OldRoot], -preserve_lengths => 1) if ($RootedTest); #Remove the old root node if the tree was rooted as it is now useless
	
	print STDERR "Rerooted tree\n" if ($RootedTest); 
	
	return(1);
}

=pod
=item * FindTrueRoot(BioPerl TreeObject, Array Pointer To Outgroup BioPerl Node IDs)
A function that will reroot a tree given an outgroup array of bioperl nodeIDs.
=cut

sub TreeIntersection($$$){
	
	my ($treeAio,$treeBio,$verboseintersection) = @_;
	 #Using IO::String to create io hadles for the newick strings. Do this externally to this function using my $io = IO::String->new($string);
	#$verboseintersection is a flag for printing out a whole load of info regarding the trees and which nodes intersect.
		
	my $input = new Bio::TreeIO(-fh   => $treeAio,
	                            -format => "newick") or die $!;
	                            
	my $TreeAObject = $input->next_tree;
	
	$input = new Bio::TreeIO(-fh   => $treeBio,
	                         -format => "newick") or die $!;
	
	my $TreeBObject = $input->next_tree;
	#Forgive the use of A and B as identifiers. I feel that this is the best way to differentiate between two very similiar objects.
	
	my @Ataxa = map{$_->id}$TreeAObject->get_leaf_nodes;
	my @Btaxa = map{$_->id}$TreeBObject->get_leaf_nodes;
	
	print STDERR "Tree A leaves: [".scalar(@Ataxa)."] \n";
	print STDERR join(',',sort(@Ataxa));
	print STDERR "\n";
	
	print STDERR "Tree B leaves: [".scalar(@Btaxa)."] \n";
	print STDERR join(',',sort(@Btaxa));
	print STDERR "\n";
	
	my ($Union,$Intersection,$ListAExclusive,$ListBExclusive) = IntUnDiff(\@Ataxa,\@Btaxa);
	
	die "No taxa in commonbetween the two trees!\n" if (scalar(@$Intersection) == 0);
	
	#A Remove
	foreach my $ANodeToRemove (@$ListAExclusive){
		
		my $NodeID = $TreeAObject->find_node(-id => $ANodeToRemove);
		$TreeAObject->remove_Node($NodeID);
	}
	
	#Tidy up tree and remove 'dangling branches'
	
	my @DanglingNodes = grep{$_->is_Leaf}$TreeAObject->find_node(-id => '');
	
	while (scalar(@DanglingNodes) > 0 ){
	
		foreach my $DanglingNode (@DanglingNodes){
	
			$TreeAObject->remove_Node($DanglingNode);
		}
			
		@DanglingNodes = grep{$_->is_Leaf}$TreeAObject->find_node(-id => '');
	}
	
	#Bremove
	foreach my $BNodeToRemove (@$ListBExclusive){
		
		my $NodeID = $TreeBObject->find_node(-id => $BNodeToRemove);
		$TreeBObject->remove_Node($NodeID);
	}
	
	#Tidy up tree and remove 'dangling branches'
	
	@DanglingNodes = grep{$_->is_Leaf}$TreeBObject->find_node(-id => '');
	
	while (scalar(@DanglingNodes) > 0 ){
	
		foreach my $DanglingNode (@DanglingNodes){
	
			$TreeBObject->remove_Node($DanglingNode);
		}
			
		@DanglingNodes = grep{$_->is_Leaf}$TreeBObject->find_node(-id => '');
	}
	
	# Output tree descriptions
	
	if($verboseintersection){ #Some useful info regaridng leaves of different trees
		print STDERR "Taxa in common: (".scalar(@$Intersection).")\n";
		my $OutString = join(',',sort(@$Intersection));
		print STDERR $OutString."\n";
		
		print STDERR "Tree A Tree Unique leaves: (".scalar(@$ListAExclusive).")\n";
		$OutString = join(',',@$ListAExclusive);
		print STDERR $OutString."\n";
		
		print STDERR "Output A Tree leaves: (".scalar($TreeAObject->get_leaf_nodes).")\n";
		$OutString = join(',',sort(map{$_->id}$TreeAObject->get_leaf_nodes));
		print STDERR $OutString."\n";
		
		print STDERR "Tree B Tree Unique leaves: (".scalar(@$ListBExclusive).")\n";
		$OutString = join(',',@$ListBExclusive);
		print STDERR $OutString."\n";
		
		print STDERR "Output B Tree leaves: (".scalar($TreeBObject->get_leaf_nodes).")\n";
		$OutString = join(',',sort(map{$_->id}$TreeBObject->get_leaf_nodes));
		print STDERR $OutString."\n";
	}
	
	return($TreeAObject,$TreeBObject);	
}

=pod
=item * TreeIntersection(BioPerl TreeObjectA, BioPerl TreeObjectB, verbose flag)
A function to take two newick tree strings (using IO::String) and output two new bioperl tree objects, both of which contain the same nodes (i.e. they can only differ on topology and branch lengths, not on leaves).

=cut

sub TreeHash2Newick($$){

	my ($TreeCacheHash,$DesiredRoot) = @_;

	print "Root: ".$DesiredRoot."\n";

	my $CladeNodes = $TreeCacheHash->{$DesiredRoot}{'all_Descendents'};
	my @rootDescendents = @{$TreeCacheHash->{$DesiredRoot}{'each_Descendent'}};;
	push(@$CladeNodes,$DesiredRoot) unless(scalar(@rootDescendents) == 1);

	my $CladeNodesHash = {};
	map{$CladeNodesHash->{$_}=undef}@$CladeNodes; # Initialise a hash with keys equal to the clade nodes. This allows for easier lookup and deleting of elements

	foreach my $CladeNode (keys(%$CladeNodesHash)){
		
		if ($TreeCacheHash->{$CladeNode}{'is_Leaf'}){
		
			my ($nodename,$edge_length) = (($TreeCacheHash->{$CladeNode}{'node_id'}),($TreeCacheHash->{$CladeNode}{'branch_length'}));
				
			my $NewickString = $nodename.":".$edge_length;
			
			$TreeCacheHash->{$CladeNode}{'Newick_clade'} = $NewickString; #Update old newick string from 'NULL' to correct value
		
			delete($CladeNodesHash->{$CladeNode});
		}
	}
	#Initilise leaf nodes with values
	
	
	#(possibility of additional functionality - exclude a list of genomes)
		
	while (scalar(keys(%$CladeNodesHash))){

		foreach my $node (keys(%$CladeNodesHash)) {
					
			my @Descendants = @{$TreeCacheHash->{$node}{'each_Descendent'}};
			die "Non Binary Tree for node: $node  id: $TreeCacheHash->{$node}{'node_id'}     - Does not have two descendents\n" unless(scalar(@Descendants) == 2 || $node eq $DesiredRoot); #Check for a binary tree
			#There's no reason that we can't have a non-binary tree, I just need to change the next few lines into loops
			
			my ($Descendent1NewickClade,$Descendent2NewickClade) = ('NULL','NULL');
				
			{no warnings 'uninitialized';
			$Descendent1NewickClade = $TreeCacheHash->{$Descendants[0]}{'Newick_clade'} if(exists($TreeCacheHash->{$Descendants[0]}{'Newick_clade'}));
			$Descendent2NewickClade = $TreeCacheHash->{$Descendants[1]}{'Newick_clade'} if(exists($TreeCacheHash->{$Descendants[1]}{'Newick_clade'}));
			#Extract the newick format trees of the children of node
			}
									
			if($Descendent1NewickClade ne 'NULL' && $Descendent2NewickClade ne 'NULL'){
				
				#Aggregate two nodes below this one into a newick format
				my ($nodename,$edge_length) = ($TreeCacheHash->{$node}{'node_id'},$TreeCacheHash->{$node}{'branch_length'});
				
				my $NewickString = "(".$Descendent1NewickClade.",".$Descendent2NewickClade.")".$nodename.":".$edge_length;
				$TreeCacheHash->{$node}{'Newick_clade'} = $NewickString; #Update old newick string from 'NULL' to correct value
				
				delete($CladeNodesHash->{$node});
			}	#i.e. if the nodes below the one under study have newick strings of the clade beneath them, then we can proceed. Else wait until the next loop
		}
	}
	

	$TreeCacheHash->{$DesiredRoot}{'Newick_clade'} = $TreeCacheHash->{$rootDescendents[0]}{'Newick_clade'} if(scalar(@rootDescendents) == 1);	
	#In some cases, the root of the tree only has one descendent (or this is just how bioperl does things). This is a quick counter
	
	##Ugly!
	my $FullTree = $TreeCacheHash->{$DesiredRoot}{'Newick_clade'};
	$FullTree="(".$FullTree.");";
	$TreeCacheHash->{$DesiredRoot}{'Newick_clade'}=$FullTree;
	#Above is to correct the final trees created in the resulting hashes.
	
	return($FullTree);
}

=pod
=item *TreeHash2Newick($rootleft,$rootright)

Given the left_id and the right_id of a node, this will return the tree stored in SUPERFAMILY under that node in the tree table as a newick string.
Both right and left ids need to be specified so as to prevent the all to easy error choosing the wrong node as the root. Returns two hashes: one has
keys of left_ids and which hash to: $LeftHash->{leftid}=[rightid,[($SQLnodename,$edge_length,$CladeAsNewickString)]]. The other hash is the same but with 
rightids.

This function only deals with binary trees, as in SUPERFAMIY.
=cut


sub FindMRCA($$$){
	
	my ($TreeCacheHash,$root,$LeavesArrayRef) = @_;
	
	my $Ancestor = $$LeavesArrayRef[0];
	my @Clade = ($Ancestor); #Initialise the clade under study with a random leaf 
	
	my $Descendent;
	
	my ($Union,$Intersection,$ListAExclusive,$ListBExclusive) = IntUnDiff($LeavesArrayRef,\@Clade);

	while (scalar(@$ListAExclusive)){
		
		die "root doesn't have all leaves submitted as descendents - error!\n" if ($Ancestor eq $root);
		
		$Descendent = $Ancestor;
		$Ancestor = $TreeCacheHash->{$Descendent}{'ancestor'};
		
		@Clade = @{$TreeCacheHash->{$Ancestor}{'all_Descendents'}};
		
		($Union,$Intersection,$ListAExclusive,$ListBExclusive) = IntUnDiff($LeavesArrayRef,\@Clade);
	}#Climb up the tree while not all leaves are embodied by the MRCA under consideration
	
	return($Ancestor);
}

=pod
=item *FindMRCA($TreeCacheHash,$root,$LeavesArrayRef)

Finds the most recent common ancestor of the leaves provided in $LeavesArrayRef. Works by climbing up the tree incrementally

=cut


1;

