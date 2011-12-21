package Supfam::TreeFuncsNonBP;
require Exporter;

=head1 NAME

Supfam::TreeFuncsNonBP

=head1 SYNOPSIS

Holds functions related to parsing Newick trees and calculating domain architecture information 
This is a branch from the TreeFuncs.pm but with the intention of it being written in pure perl - i.e. no bioperl! It doesn't really add all that much except cut down on dependencies and increase code transparency

use Supfam::TreeFuncsNonBP;

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
			FindMRCA
			Newick2Node
			DolloPLeavesWithTrait
			RAxMLLeavesWithTrait
			Node2Newick
			Splice_Node
			Splice_Root
			GraphPartition
			Node2Newick
			ExtractNewickSubtree
			RootByOutgroup
			RootOnInternalNode
			AllAncestors
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use DBI;

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
	
	my ($NewickTreeString) = @_;
	
	my $tic = Time::HiRes::time; 
	
	my $TreeHash = {}; # A massive limitation of this script is the use of BioPerl TreeIO which is SLOOOOW. This is a lookup hash to speed things up.
	die "No semi-colon on tail of Newick string $! $?\n" unless($NewickTreeString =~ m/;$/); #A little error checking
	chomp($NewickTreeString);
	
	
	#Read in and initialise tree
	
	my $root = Newick2Node($NewickTreeString,$TreeHash);
	#Recursive function that populates tree entries into $TreeCacheHash
	
	my $BranchLengthsFlag = 1; #Flag to test if the input file contained branch length info (will no compute cladetime hash if not)
	
	foreach my $node (keys(%$TreeHash)){
		
		if ($TreeHash->{$node}{'branch_length'} ~~ undef){$BranchLengthsFlag =0; }
		my @CladeLeavesList;
		
		$TreeHash->{$node}{'Clade_Leaves'} = \@CladeLeavesList;
		
		@CladeLeavesList = grep{$TreeHash->{$_}{'is_Leaf'} == 1}(@{$TreeHash->{$node}{'all_Descendents'}});	
	}
	
	# Add $TreeHash->{$node}{'Clade_Leaves'}
	
#	EasyDump("./Treedraft",$TreeHash);
	
	if($BranchLengthsFlag){
		
		map{GenerateCladeTimeHash($_,$TreeHash)}(@{$TreeHash->{$root}{'all_Descendents'}},$root);
	}else{
		
		$TreeHash->{$root}{'branch_length'} = undef if($TreeHash->{$root}{'branch_length'}); #This is a hideous way to deal with the fact that it is possible for the root to have a branch length when no other node does.
		
		print STDERR "Branch lengths not specieifed in input tree\n";
	}
	#This adds two other key/value pairs:
	# 	$TreeHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
	# 	$TreeHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
	
	#$BranchLengthsFlag is only for the case where the input tree has no branch lengths;
	
	my $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Non-BioPerl Tree Cache hash:".($toc-$tic)."seconds\n";
	
	return($root,$TreeHash);
}

=pod
=item * BuildTreeCacheHash($NewickTreeString)
Often it is much simpler to simply provide a string of a newick tree to this subroutine and recieve a hash of that tree with a whole load
of information about each node precalculated.

Each node is stored as follows $TreeHash->{Node key}
There is then a variety of values stored:
$TreeHash->{NodeID}{'all_Descendents'} - All Node keys of nodes below this one, stored as a hash {node}=undef to allow for quick querying
$TreeHash->{NodeID}{'each_Descendent'} - The direct descendents of this node
$TreeHash->{NodeID}{'Clade_Leaves'} - All the leaf nodes below this one, stored as a hash {leaf}=undef to allow for quick querying
$TreeHash->{NodeID}{'branch_length'} - The lengDolloParsimonyAncestralStateth of the branch connecting this node to it's parent
$TreeHash->{NodeID}{'is_Leaf'} - 1/0 flag for if the node is a leaf
$TreeHash->{NodeID}{'ancestor'} - BioPerl reference of nodes parent
$TreeCacheHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
$TreeCacheHash->{$MRCA}{'node_id'}='Node name' as given in the input newick file

returned is:

($root,$TreeCacheHash);

$root - the hash key of the root (remeber that the keys are the clades in newick)

$TreeCacheHash - as above

=cut


sub Newick2Node{
	
	#Provide a node and the tree hash, return a pointer to the hash of {cladepoint => node}. Cladepoint is equal to the normalised (0<x<1) point in the total clade time at which the node sits
	my ($NewickString,$TreeHash,$Ancestor) = @_;
	
	$Ancestor = 'ROOT' if (@_ < 3); #If $Ancestor isn't set, we shall assume that it's the root of the tree


##Construct regexes that parse the newick format - these took a day to make!! I could move these outside of this sub, but speed really isn't an issue.
our $CladeREGEX;
#Example string: ^A:0.1,B:0.2,(C:0.3,D:0.4):0.5$     or      (zf:0.038180391033334815,(ML:0.03567456015116893,gg:0.02024961624308485):0.008713385399205688):0.33501485928240105

$CladeREGEX = qr{ 
(?:
		(
		      (?: \w+ (?: :\d*\.?\d*(?:[eE][-+][0-9]+)?)?   )
		  #Match simple clades of one item like: A:0.1 -  B:0.2
		| 
		     (:? \( (??{$CladeREGEX}) \) \w* (?: :\d*\.?\d*(?:[eE][-+][0-9]+)?)?   ) 
		 #Initialise $CladeREGEX again so as to parse more complex clades, e.g (C:0.3,D:0.4):0.5
		)
			(?:  , (??{$CladeREGEX})  )?
	)
}x;
	

our $SubcladeRegex;
#Have to predeclare else variable wont be in scope (it calls itself)
$SubcladeRegex= qr{ #Example string: ^A:0.1,B:0.2,(C:0.3,D:0.4):0.5$     or      (zf:0.038180391033334815,(ML:0.03567456015116893,gg:0.02024961624308485):0.008713385399205688):0.33501485928240105
		\( ((??{$CladeREGEX})) \) (\w*) (?: :(\d*\.?\d* (?:[eE][-+][0-9]+)?  ))?  #Initialise $CladeREGEX again so as to parse more complex clades, e.g (C:0.3,D:0.4):0.5
}x;


our $UngroupedSubcladeRegex;#Have to predeclare else variable wont be in scope (it calls itself)
$UngroupedSubcladeRegex= qr{
	(?:
		(
			   \w+ (?: :\d*\.?\d* (?:[eE][-+][0-9]+)?  )?
		|
			  \( (??{$CladeREGEX}) \) \w* (?: :\d*\.?\d* (?:[eE][-+][0-9]+)?  )?
		)
		,?
	)
 }x;



	
if($NewickString =~ m/;$/){ #Deal with full newick strings

	if($NewickString =~ m/^\(?($SubcladeRegex)\)?;$/){
		
		$NewickString = $1;
		
		Newick2Node($NewickString, $TreeHash,'ROOT');

		return($NewickString); #i.e. returns the root oDolloParsimonyAncestralStatef the tree node in the hash
		
	}elsif($NewickString =~ m/^\(?($CladeREGEX)\)?;$/){
		
		#Another popular style to provide a tree with no distance to a root: (A,B,(C,D)); :0.1,:0.2,(:0.3,:0.4):0.5);  - we shall trim off the semi-colon and parse through (:0.1,:0.2,(:0.3,:0.4):0.5):0.0
		
		$NewickString = "(".$1."):0.00";
		#this could end up adding a zero length branch length to the root in a tree without branch lengths. At this point, we don't know if the tree has branch lenghts or not, but this will be checked for in the parent 'BuildTreeCacheHash' routine
		
		Newick2Node($NewickString, $TreeHash,'ROOT');	
		
		return($NewickString); #i.e. returns the root of the tree node in the hash
	
	}else{
		
		die "Unable to parse tree!!!! \n";
	}

}

my $All_Descendants_Array_Ref = []; #Decendants will be stored as an array ref

	if($NewickString =~ m/^$SubcladeRegex$/){
		# Node name optional, but a branch length is specified, example: (:0.1,:0.2,(:0.3,:0.4):0.5):0.0  or (A:0.1,B:0.2,(C:0.3,D:0.4):0.5)ROOT:0.0
		# OR node name given, but a branch length is not specified, example: (A,B,(C,D)) or (A,B,(C,D)E)F (i.e only leaf nodes labelled and all nodes labelled)

		my ($DescendentString,$NodeName,$Branchlength) = ($1,$2,$3);
		
		#print $NewickString."  = Full \n";
		#print $DescendentString;
		#print "   - String \n";

		$TreeHash->{$NewickString}{'is_Leaf'} = 0;
		$TreeHash->{$NewickString}{'branch_length'} = $Branchlength;
		$TreeHash->{$NewickString}{'node_id'} = $NodeName;
		
		$TreeHash->{$NewickString}{'ancestor'} = $Ancestor;
		
		#$DescendentString looks like A:0.1,B:0.2,(C:0.3,D:0.4):0.5  - so splitting on comma will give each desendent clade
		
		my @DirectDescendentClades;
		
		while($DescendentString =~ m/$UngroupedSubcladeRegex/g){push(@DirectDescendentClades,$1)};

		#print join("\n",@DirectDescendentClades);
		#print "   - Clade \n";

	##might kickDolloParsimonyAncestralState
	
		$TreeHash->{$NewickString}{'each_Descendent'} = \@DirectDescendentClades;
		$TreeHash->{$NewickString}{'all_Descendents'} = $All_Descendants_Array_Ref;
		
		@$All_Descendants_Array_Ref = @DirectDescendentClades;
		
		foreach my $DescendentClade (@DirectDescendentClades){
		
			my $ChildCladeAllDescendantsArrayRef = Newick2Node($DescendentClade,$TreeHash,$NewickString);
			push(@$All_Descendants_Array_Ref,@$ChildCladeAllDescendantsArrayRef);
		}
		
	}elsif($NewickString =~ m/^(\w+):?(\d+\.?\d+(?:[eE][-+][0-9]+)?)?$/){
		
		#Clade is in fact a leaf like A:0.1 or :0.1 or A
		
		my ($NodeName,$Branchlength) = ($1,$2);
		
		#print "$NewickString\n";
		
		$TreeHash->{$NewickString}{'is_Leaf'} = 1;
		$TreeHash->{$NewickString}{'branch_length'} = $Branchlength;
		$TreeHash->{$NewickString}{'node_id'} = $NodeName;
		$TreeHash->{$NewickString}{'each_Descendent'} = [];
		$TreeHash->{$NewickString}{'all_Descendents'} = $All_Descendants_Array_Ref;
		$TreeHash->{$NewickString}{'ancestor'} = $Ancestor;
		
	}else{
		
		die "Cannot parse tree string. Perhaps the tree is in an unusual format? $NewickString \n";	
	}
	
	return($All_Descendants_Array_Ref);
}

=pod
=item * Newick2Node($MRCA, $TreeCacheHash, [Ancestor])
An internal function for BuildTreeCacheHash. Adds these additional values onto hash:

#Each node is stored as follows $TreeHash->{CladeID}
#$TreeHash->{CladeID}{'all_Descendents'} -DolloParsimonyAncestralState All hash keys of nodes below this one
#$TreeHash->{CladeID}{'each_Descendent'} - The direct descendents of this node (again an array ref to a list of hash keys)
#$TreeHash->{CladeID}{'is_Leaf'} - 1/0 flag for if the node is a leaf
#$TreeHash->{CladeID}{'ancestor'} - hash key of node parent
#$TreeCacheHash->{CladeID}{'node_id'}='Node name' as given in the input newick file
#$TreeCacheHash->{CladeID}{'branch_length'} = length of ancestral branch to parent node

Note that if these values weren't specified in the newick input file, these values will be undef

Don't call this function directly please.
=cut


sub GenerateCladeTimeHash($$){
	
	#Provide a node and the tree hash, return a pointer to the hash of {cladepoint => node}. Cladepoint is equal to the normalised (0<x<1) point in the total clade time at which the node sits
	my ($MRCA, $TreeCacheHash) = @_;
	#MRCA is the most recent common acnestor of clade under study

	#Total Branch Length in clade below this point
	my $TotalBranch = 0;
	map{$TotalBranch += $TreeCacheHash->{$_}{'branch_length'}}(@{$TreeCacheHash->{$MRCA}{'all_Descendents'}}); #Sum all branch lengths beneath MRCA of clade
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
An internal function for BuildTreeCacheHash. Adds two additiprint STDERR "WARNING: Branch Lengths desired, but none givne in input tree\n" unless($TreeHash->{$desired_root}{'branch_length'});onal values onto hash:
$TreeCacheHash->{CladeID}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
Don't call this function directly please.
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

sub ExtractNewickSubtree($$$$){
	
	my ($TreeHash, $desired_root,$branchesflag,$internalnodesflag) = @_;

	die "Desired root not in the tree!" unless(exists($TreeHash->{$desired_root})); 
	#Error checking
	
	my @Descendents = @{$TreeHash->{$desired_root}{'each_Descendent'}};

	my @Subclades;
	
	foreach my $Descendent (@Descendents){
		
		my $SubCladeInNewick = Node2Newick($TreeHash,$Descendent,$desired_root,$branchesflag,$internalnodesflag);
		push(@Subclades,$SubCladeInNewick);
	}
	
	my $CladeInNewick = "(".join(',',@Subclades).")";
		
	if($internalnodesflag){
	
			print STDERR "WARNING: internal node_ids desired, but none givne in input tree\n" unless($TreeHash->{$desired_root}{'node_id'});
			my $NodeName = $TreeHash->{$desired_root}{'node_id'};
			$CladeInNewick = $CladeInNewick.$NodeName;
	}
		
	if($branchesflag){
	
			print STDERR "WARNING: Branch Lengths desired, but none givne in input tree\n" unless($TreeHash->{$desired_root}{'branch_length'});
			my $BranchLength = $TreeHash->{$desired_root}{'branch_length'};
			$CladeInNewick = $CladeInNewick.':'.$BranchLength;
	}
	
	$CladeInNewick = $CladeInNewick.";";
	
	return($CladeInNewick);
}

=pod
=item *ExtractNewickSubtree($TreeHash, $desired_root,$branchesflag,$internalnodesflag)

Given a tree as a tree hash $TreeHash and a desired root, this function will give the tree beneath that node in newick. 
You can specify how you want the tree using $branchesflag (1 - include branch lengths) and $internalnodesflag (1 - include the names of internal nodes).  Leaf names are included automatically

=cut

sub RootOnInternalNode{
	
	my ($TreeHash, $desired_root,$branchesflag,$internalnodesflag) = @_;
	
	die 'Desired root is already the root! Dont make calls to useless subroutines when an if($node == $root) statement is sufficient!'."\n\n" if($TreeHash->{$desired_root}{'ancestor'} eq 'ROOT');
	die "Desired root not in the tree!" unless(exists($TreeHash->{$desired_root})); 
	#Error checking
	
	my @Descendents = @{$TreeHash->{$desired_root}{'each_Descendent'}};
	
	my @NeighbouringNodes = (@Descendents,$TreeHash->{$desired_root}{'ancestor'});
	
	my @Subclades;
	
	foreach my $Neighbour (@NeighbouringNodes){
		
		my $SubCladeInNewick = Node2Newick($TreeHash,$Neighbour,$desired_root,$branchesflag,$internalnodesflag);
		push(@Subclades,$SubCladeInNewick);
	}
	
	my $CladeInNewick = "(".join(',',@Subclades).")";
		
	if($internalnodesflag){
	
			print STDERR "WARNING: internal node_ids desired, but none givne in input tree\n" unless($TreeHash->{$desired_root}{'node_id'});
			my $NodeName = $TreeHash->{$desired_root}{'node_id'};
			$CladeInNewick = $CladeInNewick.$NodeName;
	}
		
	$CladeInNewick = $CladeInNewick.";";
	
	return($CladeInNewick);	
}

=pod
=item *RootOnInternalNode($TreeHash,$desired_root,$branchesflag,$internalnodesflag)

Given a tree as a tree hash $TreeHash and a desired root, this function will give that tree rooted on the desired node in newick. 
You can specify how you want the tree using $branchesflag (1 - include branch lengths) and $internalnodesflag (1 - include the names of internal nodes). Leaf names are included automatically

=cut

sub RootByOutgroup{
	
	my ($TreeHash,$Root,$OutgroupLeaves,$branchesflag,$internalnodesflag,$RootName,$RootBranchLength) = @_ ;
	
	my ($RootAncestor,$RootDescendant) = FindTrueRoot($TreeHash,$OutgroupLeaves,$Root);
	
	my $CladeInNewick;
	
	unless($RootAncestor eq $Root){
	
		Splice_Root($TreeHash,$Root);
		
		my @Subclades;

		my $SubCladeInNewick = Node2Newick($TreeHash,$RootAncestor,$RootDescendant,$branchesflag,$internalnodesflag);
		push(@Subclades,$SubCladeInNewick);
	
		$SubCladeInNewick = Node2Newick($TreeHash,$RootDescendant,$RootAncestor,$branchesflag,$internalnodesflag);
		push(@Subclades,$SubCladeInNewick);
		
		$CladeInNewick = "(".join(',',@Subclades).")";
			
		if($internalnodesflag){

				my $NodeName = $RootName;
				$CladeInNewick = $CladeInNewick.$NodeName;
		}
			
		if($branchesflag){
		
				my $BranchLength = '0:00';
				$CladeInNewick = $CladeInNewick.':'.$BranchLength;
		}
		
		$CladeInNewick = $CladeInNewick.";";
		
	}else{
		
		print STDERR "Tree is already rooted on the appropriate root\n";
		
		$TreeHash->{$Root}{'branch_length'} = $RootBranchLength if($branchesflag);
		$TreeHash->{$Root}{'node_id'} = $RootName if($internalnodesflag);
		
		$CladeInNewick = ExtractNewickSubtree($TreeHash, $RootAncestor,$branchesflag,$internalnodesflag);
	}
	
	undef $TreeHash; #The hash is no good to us after we have spliced the root out of it!
	
	return($CladeInNewick);
}

=pod
=item *RootAtMidpoint($TreeHash,$branchesflag,$internalnodesflag)

Given a tree as a tree hash $TreeHash and a desired root, this function will give that tree rooted on the desired node in newick. 
You can specify how you want the tree using $branchesflag (1 - include branch lengths) and $internalnodesflag (1 - include the names of internal nodes). Leaf names are included automatically

=cut

sub Splice_Node($$){
	
	my ($TreeHash,$Node2Splice) = @_;
	
	my $Ancestor = $TreeHash->{$Node2Splice}{'ancestor'};
	my @Descendants = @{$TreeHash->{$Node2Splice}{'each_Descendent'}};
	
	die "Can't Splice_midpoint on root - use Splice_Root for that\n" if($Ancestor eq 'ROOT');

	@{$TreeHash->{$Ancestor}{'each_Descendent'}} = grep{!/$Node2Splice/}@{$TreeHash->{$Ancestor}{'each_Descendent'}};
	#Remove $Node2Splice from ancesotrs 'each_Descendant array record

	foreach my $SingleDescendant (@Descendants){
	
		my $TotalBranchLength = ($TreeHash->{$SingleDescendant}{'branch_length'}) + ($TreeHash->{$Node2Splice}{'branch_length'});

		$TreeHash->{$SingleDescendant}{'branch_length'} = $TotalBranchLength;
		$TreeHash->{$SingleDescendant}{'ancestor'} = $Ancestor; #i.e. change the ancestor of the descendant to the parent node.
	
		push(@{$TreeHash->{$Ancestor}{'each_Descendent'}},$SingleDescendant); #Add single_descendant onto the ancesotrs direct descendants list
	}
	
	
	#Remove $Node2Splice from all of it's ancestors 'all_descendants' arrays
	my $IncrementalAncestor = $Ancestor;
	
	while ($IncrementalAncestor ne 'ROOT'){
		
		my $AllDescendants = $TreeHash->{$IncrementalAncestor}{'all_Descendents'};
		@$AllDescendants = grep{!/$Node2Splice/}$AllDescendants;
		
		@{$TreeHash->{$IncrementalAncestor}{'Clade_Leaves'}} = grep{!/$Node2Splice/}@{$TreeHash->{$IncrementalAncestor}{'Clade_Leaves'}} if($TreeHash->{$IncrementalAncestor}{'is_Leaf'});
		#Removes Leaf entries if the node being spliced is a leaf
		
		$IncrementalAncestor = $TreeHash->{$IncrementalAncestor}{'ancestor'};
	}
	
	
			
	delete $TreeHash->{$Node2Splice}; #Finally, delete the spliced hash node
	
	return(1);
}

=pod
=item *Splice_Node($TreeHash,$Node2Splice)

Removes a non-root node (regardless of the degree of the node). Branch lengths are updated accordingly, as are 'ancestor', 'each_Descendant' and 'all_Decendants' entries.

=cut

sub Splice_Root($$){
	
	my ($TreeHash,$Root2Splice) = @_;
	
	die "This is not the root of the tree! Use Splice_Midpoint sub for removing a midpoint node of degree 2" if($TreeHash->{$Root2Splice}{'ancestor'} eq 'ROOT');
	
	my @Descendants = (@{$TreeHash->{$Root2Splice}{'each_Descendent'}});
	
	die "This sub is for removing root nodes with two neighbours ONLY! If your tree has three (or more) descendants, then it is rooted on an internal node, so don't remove it\n" unless(scalar(@Descendants) == 2 || $TreeHash->{$Root2Splice}{'is_Leaf'});
	
	my ($DescendantA, $DescendantB) = @Descendants;
	
	my $TotalBranchLength = ($TreeHash->{$DescendantA}{'branch_length'}) + ($TreeHash->{$DescendantB}{'branch_length'});
	
	$TreeHash->{$DescendantA}{'ancestor'} = $DescendantB;
	$TreeHash->{$DescendantB}{'ancestor'} = $DescendantA;
	
	$TreeHash->{$DescendantA}{'branch_length'} = $TotalBranchLength;
	$TreeHash->{$DescendantB}{'branch_length'} = $TotalBranchLength;
	
	delete $TreeHash->{$Root2Splice}; #Finally, delete the spliced root
	
	#NOTE We now effectively have two seperate trees backing on to each other. We have a graph! Not a tree, as we have no root! Stop thinking about this as a tree until you re root it
	# I wouldn't call any other tree subs on this hash. Reroot it (which will give you a fresh newick string) and then create a fresh tree hash.  In fact, don't use this sub unless you really understand it
	
	return(1);	
}

=pod
=item *Splice_Root($TreeHash,$Root2Splice)

If a tree is rooted at a midpoint between two nodes, this will remove that node and join the two descendents (i.e. they each will consider each other the others ancestor, with branch lengths udated accordingly).
This breaks the idea of a tree; it's now a graph. I wouldn't call any tree subs on this hash - reroot it and create a new hash (deleting this hash so as to save on memory).
 
This is a semi-internal sub. Don't call it directly unless you fully understand what it does.

=cut


sub Node2Newick($$$$$);

sub Node2Newick($$$$$){
	
	my ($TreeHash,$Ancestor,$node,$branchesflag,$internalnodesflag) = @_;
	
	my $CladeInNewick;
	my $BranchLength;	

	if($branchesflag){
	
		print STDERR "WARNING: Branch Lengths desired, but none givne in input tree\n" unless($TreeHash->{$node}{'branch_length'});
			
		unless(grep{/$Ancestor/}(@{$TreeHash->{$node}{'each_Descendent'}})){ #Unless the desired ancestor is currently a descendant of the node
				
			$BranchLength = $TreeHash->{$node}{'branch_length'};
				
		}else{
				
			$BranchLength = $TreeHash->{$Ancestor}{'branch_length'}; #If the desired ancestor is a descendant, use the branch between the two
		}
			
	}	
	#Correctly sets branch lengths if they were requested. Can deal with the desired ancestor currently being a descendant
	
	
	unless($TreeHash->{$node}{'is_Leaf'}){
		
		my @NeighbouringNodes = (@{$TreeHash->{$node}{'each_Descendent'}},$TreeHash->{$node}{'ancestor'}); #Forget about direction, we're treating the tree as a graph (no directionaility yet)
		my @NewDescendants = grep{$_ ne $Ancestor}@NeighbouringNodes; #Specifying an ancestory gives directionality to the tree - non-ancestor neighbours are therefore descendents
		
		print scalar(@NeighbouringNodes)."\n";
		
		my @SubClades;
		
		foreach my $Descendent (@NewDescendants){
			
			my $SubCladeInNewick = Node2Newick($TreeHash,$node,$Descendent,$branchesflag,$internalnodesflag);
			push(@SubClades,$SubCladeInNewick);
		}
		
		$CladeInNewick = "(".join(',',@SubClades).")";
		
		$CladeInNewick = $CladeInNewick.':'.$BranchLength if($branchesflag);	
		
	}else{ #i.e. it's a leaf
		
		my $NodeName = $TreeHash->{$node}{'node_id'};
		die "WARNING: no node_id for leaf in input tree\n" unless($TreeHash->{$node}{'node_id'});
		
		$CladeInNewick = $NodeName;
		
		$CladeInNewick = $CladeInNewick.':'.$BranchLength if($branchesflag);	
		
	}
	
	return($CladeInNewick);	
}

=pod
=item *Node2Newick($TreeHash,$Ancestor,$node,$branchesflag,$internalnodesflag)

Another recursive function - this time to generate the tree in newick beneath this node. You can specify how you want the tree using $branchesflag (1 - include branch lengths) and $internalnodesflag (1 - include the names of internal nodes).
Note that this function alls you to specify and ancestor of a node, allowing you to change the direction of the tree (useful when re rooting).

=cut

sub AllAncestors($$){
	
	my ($TreeHash, $Node) = @_;
	
	my $IncrementalAncestor = $TreeHash->{$Node}{'ancestor'};
	
	my $AllAncestors = [];
		
	while ($IncrementalAncestor ne 'ROOT'){
	
		push(@$AllAncestors,$IncrementalAncestor);
		
		$IncrementalAncestor = $TreeHash->{$IncrementalAncestor}{'ancestor'};
	}	
	
	return($AllAncestors);
}

=pod
=item *AllAncestors($TreeHash, $Node)

Given a tree as a treehash and a node, this function will provide a list of all ancestors of a node. Might consider making a Treehash entry for this sometime.

=cut

sub FindTrueRoot($$$) {
	
	my ($treehash,$OutgroupLeaves,$root) = @_;
	
	die "Root is on a leaf within th outgroup. What the hell kinda tree are you dealing with!?\n" if(grep{/$root/}@$OutgroupLeaves);
	
	my @All_Leaves = @{$treehash->{$root}{'Clade_Leaves'}};
	
	my (undef,undef,$IngroupLeaves,undef) = IntUnDiff(\@All_Leaves,$OutgroupLeaves);
	
	my $CurrentGenNode = $$OutgroupLeaves[POSIX::ceil(rand(@$OutgroupLeaves))];#Set as a random leaf
	my $PreviousGenNode = 'NULL-NODE-STRING';
	
	my $ExitFlag = 0;
	
	while(! $ExitFlag){
		
		my @Neighbours = (@{$treehash->{$CurrentGenNode}{'each_Descendent'}},$treehash->{$CurrentGenNode}{'ancestor'});
		@Neighbours = grep{!/$PreviousGenNode/}@Neighbours; #Stops us entering an infinite loop
			
		my $NextGenNode;
		my $BestCount=-1;
		
		foreach my $Neighbour (@Neighbours){
			
			my $PartitionNodes =[];	
			GraphPartition($treehash,$CurrentGenNode,$Neighbour,$PartitionNodes);
		
			my (undef,$IngroupIntersection,undef,undef) = IntUnDiff(\@All_Leaves,$OutgroupLeaves);
			#Intersection is the list of leaves that are in the current partition and in the ingroup
			
			my (undef,$OutgroupIntersection,undef,undef) = IntUnDiff(\@All_Leaves,$OutgroupLeaves);
			
			$NextGenNode = $Neighbour if($BestCount < scalar(@$IngroupIntersection));
			die "A node in the tree has the same number of ingroup nodes either side of it - poorly chosen outgroup!\n" if($BestCount == scalar(@$IngroupIntersection)); #Catch badly chosen outgroups
			$BestCount = scalar(@$IngroupIntersection) if($BestCount < scalar(@$IngroupIntersection));
			#We are looking for the partition that has the largest number of ingroup nodes branched away from the outgroup nodes
	
			unless(scalar(@$OutgroupIntersection)){ #i.e. if the current partition has no members of the outgroup, then we have found the new root!
				
				$ExitFlag=1;
				last; #Leave the for loop
			}
		}
		
		$PreviousGenNode = $CurrentGenNode;
		$CurrentGenNode = $NextGenNode;
	}
	
	my ($RootAncestor,$RootDescendant)= ($CurrentGenNode,$PreviousGenNode);
	
	return($RootAncestor,$RootDescendant);
}

=pod
=item * FindTrueRoot(BioPerl TreeObject, Array Pointer To Outgroup BioPerl Node IDs)
A function that will find two nodes with which to root the current tree.
=cut


sub GraphPartition($$$$);

sub GraphPartition($$$$){
	
	my ($TreeHash,$Ancestor,$node,$SubGraphNodes) = @_;
	
	unless($TreeHash->{$node}{'is_Leaf'}){
		
		my @NeighbouringNodes = (@{$TreeHash->{$node}{'each_Descendent'}},$TreeHash->{$node}{'ancestor'}); #Forget about direction, we're treating the tree as a graph (no directionaility yet)
		my @NewDescendants = grep{!/$Ancestor/}@NeighbouringNodes; #Specifying an ancestory gives directionality to the tree - non-ancestor neighbours are therefore descendents
		
		foreach my $Descendent (@NewDescendants){
			
			GraphPartition($TreeHash,$node,$Descendent,$SubGraphNodes);
		}
		
		push(@$SubGraphNodes,@NewDescendants)
		
	}
	
	
	
	return(1);	
}

=pod
=item *GraphPartition($TreeHash,$Ancestor,$node,$SubGraphNodes)

If we consider the tree a graph, this function will tell you all nodes beneath a particular node, coming in a particular direction (specified by $Ancestor). These are pushed onto ($TreeHash,$Ancestor,$node,$SubGraphNodes)
which should be an array ref.

=cut


sub TreeIntersection($$$){
	
	my ($treeAio,$treeBio,$verboseintersection) = @_;
	 #Using IO::String to create io hadles for the newick strings. Do this externally to this function using my $io = IO::String->new($string);
	#$verboseintersection is a flag for printing out a whole load of info regarding the trees and which nodes intersect.
	
	die "Sub undergoing a rewrite so as to leave the bioperl paradigm\n";
		
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


sub FindMRCA($$$){
	
	my ($TreeCacheHash,$root,$LeavesArrayRef) = @_;
	
	my $Ancestor = $$LeavesArrayRef[0];
	my @Clade = ($Ancestor); #Initialise the clade under study with a random leaf 
	
	my $Descendent;
	
	my (undef,undef,$ListAExclusive,undef) = IntUnDiff($LeavesArrayRef,\@Clade);

	while (scalar(@$ListAExclusive)){
		
		die "root doesn't have all leaves submitted as descendents - error!\n" if ($Ancestor eq $root);
		
		$Descendent = $Ancestor;
		$Ancestor = $TreeCacheHash->{$Descendent}{'ancestor'};
		
		@Clade = @{$TreeCacheHash->{$Ancestor}{'all_Descendents'}};
		
		(undef,undef,$ListAExclusive,undef) = IntUnDiff($LeavesArrayRef,\@Clade);
	}#Climb up the tree while not all leaves are embodied by the MRCA under consideration
	
	return($Ancestor);
}

=pod
=item *FindMRCA($TreeCacheHash,$root,$LeavesArrayRef)

Finds the most recent common ancestor of the leaves provided in $LeavesArrayRef. Works by climbing up the tree incrementally

=cut

sub RemoveNode{
	
	my ($TreeCacheHash,$root,$NodeToRemove) = @_;	
	
	#Check if root - do nothing if so and return
	
	#Print a warning if other methids have been called on the tree - advise them to rerun now that the node content is different
	
	#Delete node entry
	#Replace ancestor entry for this node with child nodes of current node (don't forget to update branch lengths)
	
	#If leaf, remove itself and check if it's ancestor had any children. If only one child, then remove it as well
	
	#Recurse back up the tree to the root, deleting this node from their 'all_Descendents' list
}


sub DolloPLeavesWithTrait($$$){
	
	my ($TreeHash,$root,$Trait) = @_;
	
	my @TreeLeaves = @{$TreeHash->{$root}{'Clade_Leaves'}};
		
	my $LeavesWithTrait = [];
		
	foreach my $Leaf (@TreeLeaves){
		
		my $TraitPosition = $TreeHash->{$Leaf}{'DolloP_Trait_String_Poistions_Lookup'}{$Trait};
		
		my $LeafTraitState = substr($TreeHash->{$Leaf}{'DolloPTraitStates'},$TraitPosition,1);
		
		push(@$LeavesWithTrait,$Leaf) if($LeafTraitState);		
	}
	
	return($LeavesWithTrait); #A list of leaves with the trait of interest, provided in the arguments
	
}

=pod
=item *DolloPLeavesWithTrait($TreeHash,$root,$Trait)

Quick function to find the leaves beneath a node ($root) that posses a given trait ($Trait). Returns an array ref of the list of nodes with the trait
=cut


sub RAxMLLeavesWithTrait($$$){
	
	my ($TreeCacheHash,$root,$Trait) = @_;
	
	my @TreeLeaves = @{$TreeCacheHash->{$root}{'Clade_Leaves'}};
		
	my $LeavesWithTrait = [];
	
	my $RAxML_Positions_Lookup = $TreeCacheHash->{$root}{'RAxML_Trait_String_Poistions_Lookup'};
	
	foreach my $Leaf (@TreeLeaves){
		
		my ($start,$end) = @{$RAxML_Positions_Lookup->{$Trait}};
			
		my $RAxMLMarginalProbabilitiesConcatenatedStringState1 = $TreeCacheHash->{$Leaf}{'RAxML_AncestralProbabilities'}[1];
		
		my $LeafTraitState = substr($RAxMLMarginalProbabilitiesConcatenatedStringState1,$start,($end-$start));
		
		
		
		push(@$LeavesWithTrait,$Leaf) if($LeafTraitState > 0.5);
	}
	
	return($LeavesWithTrait); #A list of leaves with the trait of interest, provided in the arguments
	
}

=pod
=item *RAxMLLeavesWithTrait($TreeHash,$root,$Trait)

Quick function to find the leaves beneath a node ($root) that posses a given trait ($Trait). Returns an array ref of the list of nodes with the trait
To be used after adding data to TreeCacheHAsh using RAxML Marginal Probabilities Parser
=cut


1;

