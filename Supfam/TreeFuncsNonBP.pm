package Supfam::TreeFuncsNonBP;
require Exporter;

=head1 NAME

Supfam::TreeFuncsNonBP

=head1 SYNOPSIS

Holds functions related to parsing Newick trees and calculating domain architecture information 
This is a branch from the TreeFuncs.pm but with the intention of it being written in pure perl - i.e. no bioperl!

It doesn't really add all that much except cut down on dependencies and increase code transparency. And it's faster! And Julian doesn't hate it!

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
			assignLeftRightIDs2TreeHash
			isrootedbinary_TreeHash
			recursivelyassigneLeftRightids
			CalculateLineage
			ExtractSubtree
			FindTrueRoot
			BuildTreeCacheHash
			GenerateCladeTimeHash
			normailseSQLTreeHashBranchForDeletions
			supfamSQL2TreeHash
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
			sanitise_TreeHash
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use lib "$ENV{HOME}/bin/perl-libs-custom";

use DBI;
use Supfam::Utils;
use Supfam::SQLFunc;
use Time::HiRes;

=pod
=head2 Methods
=over 4
=cut

sub supfamSQL2TreeHash($$){
		
	my ($rootleft,$rootright) = @_;
	
	my $tic = Time::HiRes::time; 
	
	my $SQLTreeCacheHash = {};
	
	my $dbh = dbConnect();
	
	my $sth = $dbh->prepare("SELECT left_id,right_id FROM tree WHERE left_id >= ? AND right_id <= ?;");
	$sth->execute($rootleft,$rootright);
	#Collect left_ids of all nodes beneth desired root - these will be used as the unique keys of the treehash
	
	while (my ($left_id,$right_id) = $sth->fetchrow_array()){
		
		$SQLTreeCacheHash->{$left_id}{'right_id'}=$right_id;
		$SQLTreeCacheHash->{$left_id}{'left_id'}=$left_id;
	}
	
	$sth = $dbh->prepare("SELECT nodename,edge_length,taxon_id FROM tree WHERE tree.left_id = ?;");
	
	my @AllLeftIDs = (keys(%$SQLTreeCacheHash));
	#All left ids in clade, including desired root
	
	foreach my $Clade_Node_leftid (@AllLeftIDs){
		
		$sth->execute($Clade_Node_leftid);
		
		my ($nodename,$edge_length,$taxon_id) = $sth->fetchrow_array();
				
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'all_Descendents'}
		=
		[grep{$_ > $SQLTreeCacheHash->{$Clade_Node_leftid}{'left_id'} 
		&& 
		($SQLTreeCacheHash->{$_}{'right_id'}) < ($SQLTreeCacheHash->{$Clade_Node_leftid}{'right_id'})}(@AllLeftIDs)];
		
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'each_Descendent'}
		=
		[grep{
			($SQLTreeCacheHash->{$_}{'left_id'}) == ($SQLTreeCacheHash->{$Clade_Node_leftid}{'left_id'}+1) 
		|| 
			($SQLTreeCacheHash->{$_}{'right_id'}) == ($SQLTreeCacheHash->{$Clade_Node_leftid}{'right_id'}-1)
			}(@AllLeftIDs)
		];
		
		#Nodes where left id of current node is one more than that in @AllLeftIDs or the right id is one less.
		
		my ($Ancestor) = grep{$_ == $SQLTreeCacheHash->{$Clade_Node_leftid}{'left_id'}-1 || ($SQLTreeCacheHash->{$_}{'right_id'}) == ($SQLTreeCacheHash->{$Clade_Node_leftid}{'right_id'}+1)}(@AllLeftIDs);
		
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'ancestor'} = $Ancestor;
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'is_Leaf'} = ($SQLTreeCacheHash->{$Clade_Node_leftid}{'right_id'} == ($SQLTreeCacheHash->{$Clade_Node_leftid}{'left_id'} + 1))?1:0;;
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'branch_length'}=$edge_length;
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'taxon_id'}=$taxon_id;
		$SQLTreeCacheHash->{$Clade_Node_leftid}{'node_id'}=$nodename;
		
		$sth->finish();
	}
	
	$SQLTreeCacheHash->{$rootleft}{'ancestor'} = 'ROOT'; # Set root as having 'ROOT' as its ancestor
	
	foreach my $LEFTID (@AllLeftIDs){
		my @NodeAllDescendents = @{$SQLTreeCacheHash->{$LEFTID}{'all_Descendents'}};
		$SQLTreeCacheHash->{$LEFTID}{'Clade_Leaves'} = [grep{$SQLTreeCacheHash->{$_}{'is_Leaf'}}(@NodeAllDescendents,$LEFTID)];
	}

	dbDisconnect($dbh);
	
	my $toc = Time::HiRes::time;
	print STDERR "Time taken to build the Non-BioPerl SQL Tree Cache hash:".($toc-$tic)."seconds\n";
	
	return($SQLTreeCacheHash);
}

=pod
=item *supfamSQL2TreeHash(root left id, root right id)

Given the left_id and the right_id of a node, this will return the tree stored in SUPERFAMILY under that node in the tree table as a TreeHash, which can be written to a string if you like.
Both right and left ids need to be specified so as to prevent the all to easy error choosing the wrong node as the root. Returns a hash with similar structure to BuildTreeCacheHash.

Cruciially, the tree is a binary tree only! Hence there are left and right id entries

Each node is stored as follows $TreeHash->{Node key}
There is then a variety of values stored:
$SQLTreeCacheHash->{NodeID}{'all_Descendents'} - All Node keys of nodes below this one, stored as a hash {node}=undef to allow for quick querying
$SQLTreeCacheHash->{NodeID}{'each_Descendent'} - The direct descendents of this node
$SQLTreeCacheHash->{NodeID}{'Clade_Leaves'} - All the leaf nodes below this one, stored as a hash {leaf}=undef to allow for quick querying
$SQLTreeCacheHash->{NodeID}{'branch_length'} - The lengDolloParsimonyAncestralStateth of the branch connecting this node to it's parent
$SQLTreeCacheHash->{NodeID}{'is_Leaf'} - 1/0 flag for if the node is a leaf
$SQLTreeCacheHash->{NodeID}{'ancestor'} - BioPerl reference of nodes parent
$SQLTreeCacheHash->{NodeID}{'node_id'}='Node name' as given in the input newick file
$SQLTreeCacheHash->{NodeID}{'left_id'} = node left id
$SQLTreeCacheHash->{NodeID}{'right_id'} = node right id

$TreeCacheHash->{$node}{'Total_branch_lengths'}  - the total branch length of the clade beneath this node
$TreeCacheHash->{$MRCA}{'Probability_Hash'}={ clade point => nodeID} - clade point is the point in the clade beneath the MRCA node which this node sits
$TreeCacheHash->{$MRCA}{'node_id'}='Node name' as given in the SQL table

The Node Keys are the left ids of the tree, but it is more reliable to use the $SQLTreeCacheHash->{NodeID}{'left_id'} entry in case the keys of the hash have been mucked about with

This function only deals with binary trees, as in SUPERFAMIY.

=cut

sub normailseSQLTreeHashBranchForDeletions($){
	
	my ($SQLTreeCacheHash) = @_;
	
	my @NodeIDs = (keys(%$SQLTreeCacheHash)); #Explicit way of extracting left ids of all the 
	
	map{die "Use this subroutine on an SQL treehash, created using supfamSQL2TreeHash\n" unless (exists($SQLTreeCacheHash->{$_}{'left_id'})) }@NodeIDs;
	#Error checking
	
	my $dbh = dbConnect();
	
	my $sth = $dbh->prepare("SELECT ancestral_info.comb_deleted FROM tree JOIN ancestral_info ON tree.left_id = ancestral_info.left_id WHERE tree.left_id = ?;");
	
	foreach my $CladeNode (@NodeIDs){
		
		my $NodeLeftID = $SQLTreeCacheHash->{$CladeNode}{'left_id'};
		$sth->execute($NodeLeftID);
		
		my ($DeletionsAlongBranch) = $sth->fetchrow_array();
		
		print STDERR $NodeLeftID."  <- Left ID  ".$DeletionsAlongBranch."  <- Deletions\n";
		
		$SQLTreeCacheHash->{$CladeNode}{'branch_length'}=$DeletionsAlongBranch;
		$sth->finish();
	}
	
	dbDisconnect($dbh);
	
	return(1);
}

=pod
=item * BuildTreeCacheHash($NewickTreeString)



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
		
		$BranchLengthsFlag = 0 if ($TreeHash->{$node}{'branch_length'} ~~ undef);
		my @CladeLeavesList;
		
		$TreeHash->{$node}{'Clade_Leaves'} = \@CladeLeavesList;
		
		@CladeLeavesList = grep{$TreeHash->{$_}{'is_Leaf'} == 1}(@{$TreeHash->{$node}{'all_Descendents'}});	
	}
	
	if($BranchLengthsFlag){
							
		foreach my $internal_node (@{$TreeHash->{$root}{'all_Descendents'}},$root){
			
			$TreeHash->{$internal_node}{'branch_length'} = 0 if($internal_node eq $root);
			
			if ($TreeHash->{$internal_node}{'branch_length'} ~~ 0 && $internal_node ne $root){
				my $NodeID = $TreeHash->{$internal_node}{'node_id'};
				print STDERR "Node id:".$NodeID." - branch length of zero found, introducing a pseudocount of 0.000001\n";
				$TreeHash->{$internal_node}{'branch_length'} = 0.00000001;
			}			
		}
				
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
	
	my $CladeREGEX;
	my $SubcladeRegex;
	my $UngroupedSubcladeRegex;
	
{no warnings; #stops a warning due to the variables being reevaled in a regex recursively
		
	#	#Example string: ^A:0.1,B:0.2,(C:0.3,D:0.4):0.5$ or (zf:0.038180391033334815,(ML:0.03567456015116893,gg:0.02024961624308485):0.008713385399205688):0.33501485928240105
	
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
		
	
	
	#Have to predeclare else variable wont be in scope (it calls itself)
	$SubcladeRegex= qr{ #Example string: ^(C:0.3,D:0.4):0.5$
			\( ((??{$CladeREGEX})) \) (\w*) (?: :(\d*\.?\d* (?:[eE][-+][0-9]+)?  ))?  #Initialise $CladeREGEX again so as to parse more complex clades, e.g (C:0.3,D:0.4):0.5
	}x;
	
	
	 #Have to predeclare else variable wont be in scope (it calls itself)
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
}

if($NewickString =~ m/;$/){ #Deal with full newick strings - so close to at the root (even if we are dealign with an unrooted tree!)

	if($NewickString =~ m/^\(?($SubcladeRegex)\)?;$/){
		
		#Another popular style to provide a tree with no distance to a root: (A,B,(C,D)); :0.1,:0.2,(:0.3,:0.4):0.5);  - we shall trim off the semi-colon and parse through (:0.1,:0.2,(:0.3,:0.4):0.5):0.0
		$NewickString = $1;
		my ($DescendentString,$NodeName,$Branchlength) = ($2,$3,$4);
		
		if($Branchlength ~~ undef){
			 
			$NewickString = $NewickString.":0.00";
		}
		
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

		$TreeHash->{$NewickString}{'each_Descendent'} = \@DirectDescendentClades;
		$TreeHash->{$NewickString}{'all_Descendents'} = $All_Descendants_Array_Ref;
		
		@$All_Descendants_Array_Ref = @DirectDescendentClades;
		
		foreach my $DescendentClade (@DirectDescendentClades){
		
			my $ChildCladeAllDescendantsArrayRef = Newick2Node($DescendentClade,$TreeHash,$NewickString);
			push(@$All_Descendants_Array_Ref,@$ChildCladeAllDescendantsArrayRef);
		}
		
	}elsif($NewickString =~ m/^(\w+):?(\d+(\.\d+)?(?:[eE][-+][0-9]+)?)?$/){
		
		#Clade is in fact a leaf like A:0.1 or A
		#$NewickString =~ m/^(\w+):?(\d+\.?\d+(?:[eE][-+][0-9]+)?)?$/)
		
		my ($NodeName,$Branchlength) = ($1,$2);

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
		
		my $SubCladeInNewick = Node2Newick($TreeHash,$desired_root,$Descendent,$branchesflag,$internalnodesflag);#($TreeHash,$Ancestor,$node,$branchesflag,$internalnodesflag)
		push(@Subclades,$SubCladeInNewick);
	}
	
	my $CladeInNewick = "(".join(',',@Subclades).")";
		
	if($internalnodesflag){
	
			print STDERR "WARNING: internal node_ids desired, but none given in input tree for nodeID $desired_root\n" unless($TreeHash->{$desired_root}{'node_id'});
			my $NodeName = $TreeHash->{$desired_root}{'node_id'};
			$CladeInNewick = $CladeInNewick.$NodeName;
	}
		
	if($branchesflag){
	
			print STDERR "WARNING: Branch Lengths desired, but none given in input tree for nodeID $desired_root\n" unless(exists($TreeHash->{$desired_root}{'branch_length'}));
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

sub Splice_Node($$$){
	
	my ($TreeHash,$Node2Splice,$root) = @_;
		
	my $Ancestor = $TreeHash->{$Node2Splice}{'ancestor'};
	my @Descendants = @{$TreeHash->{$Node2Splice}{'each_Descendent'}};
	
	die "Can't Splice_Node on root - use Splice_Root for that\n" if($Ancestor eq 'ROOT');

	@{$TreeHash->{$Ancestor}{'each_Descendent'}} = grep{$_ ne $Node2Splice}@{$TreeHash->{$Ancestor}{'each_Descendent'}};
	#Remove $Node2Splice from ancesotrs 'each_Descendant array record

	foreach my $SingleDescendant (@Descendants){
	
	unless ($TreeHash->{$SingleDescendant}{'branch_length'} ~~ undef){
	
		my $TotalBranchLength = ($TreeHash->{$SingleDescendant}{'branch_length'}) + ($TreeHash->{$Node2Splice}{'branch_length'});
		$TreeHash->{$SingleDescendant}{'branch_length'} = $TotalBranchLength;
	}
	
		$TreeHash->{$SingleDescendant}{'ancestor'} = $Ancestor; #i.e. change the ancestor of the descendant to the parent node.
	
		push(@{$TreeHash->{$Ancestor}{'each_Descendent'}},$SingleDescendant); #Add single_descendant onto the ancesotrs direct descendants list
	}
	
	
	#Remove $Node2Splice from all of it's ancestors 'all_descendants' arrays
	
	foreach my $IncrementalAncestor (@{$TreeHash->{$root}{'all_Descendents'}},$root){
		
		@{$TreeHash->{$IncrementalAncestor}{'all_Descendents'}} = grep{$_ ne $Node2Splice}@{$TreeHash->{$IncrementalAncestor}{'all_Descendents'}};
		@{$TreeHash->{$IncrementalAncestor}{'Clade_Leaves'}} = grep{$_ ne $Node2Splice}@{$TreeHash->{$IncrementalAncestor}{'Clade_Leaves'}} unless($TreeHash->{$IncrementalAncestor}{'is_Leaf'});
		#Removes Leaf entries if the node being spliced is a leaf
	}
	
	delete $TreeHash->{$Node2Splice}; #Finally, delete the spliced hash node
	
	return(1);
}

=pod
=item *Splice_Node($TreeHash,$Node2Splice)

Removes a non-root node (regardless of the degree of the node). Branch lengths are updated accordingly, as are 'ancestor', 'each_Descendant' and 'all_Decendants' entries.
Can be used on tree leaves as well.

=cut

sub sanitise_TreeHash{
	
	my ($TreeHash,$root,$verbose) = @_;
	
	$verbose = 1 if ($verbose ~~ undef);
	
	
	#This will be an evolving function.Good practice will be to call this on any function that modifies a tree object
	
	my @InternalTreeNodes;
	
	my $Changes =1;
	
	while ($Changes){
	
	$Changes=0; #Changes tracks if there has been any changes to the tree, allowing us to check again and see if we need to sanitise more nodes
	@InternalTreeNodes = grep{! $TreeHash->{$_}{'is_Leaf'}}@{$TreeHash->{$root}{'all_Descendents'}};
		
		foreach my $NonRootInternalNode (@InternalTreeNodes){
			
			my $NumberOfDescendents = scalar(@{$TreeHash->{$NonRootInternalNode}{'each_Descendent'}});
			Splice_Node($TreeHash,$NonRootInternalNode,$root) if($NumberOfDescendents == 1 || $NumberOfDescendents == 0); #i.e if the node one has one ancestor (is a pointless node), remove it
			
			if($verbose){
				print STDERR "Spliced $NonRootInternalNode \n" if($NumberOfDescendents == 1); 
				print STDERR "Removed $NonRootInternalNode \n" if($NumberOfDescendents == 0); 
			}
			
			$Changes++ if($NumberOfDescendents == 1 || $NumberOfDescendents == 0);
		}
	}
	
	return(1);
}

=pod
=item *sanitise_TreeHash($TreeHash,$root,$verboseflag)

Tidies up tree, removing (splicing) any non-leaf nodes with only a single descendent.

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
	
		print STDERR "WARNING: Branch Lengths desired, but none givne in input tree\n" unless(exists($TreeHash->{$node}{'branch_length'}));
			
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

		my @SubClades;
		
		foreach my $Descendent (@NewDescendants){
			
			my $SubCladeInNewick = Node2Newick($TreeHash,$node,$Descendent,$branchesflag,$internalnodesflag);
			push(@SubClades,$SubCladeInNewick);
		}
		
		$CladeInNewick = "(".join(',',@SubClades).")";
		$CladeInNewick = $CladeInNewick.':'.$BranchLength if($branchesflag);	
		
	}else{ #i.e. it's a leaf
		
		my $NodeName = $TreeHash->{$node}{'node_id'};
		die "WARNING: no node_id for leaf in input tree\n" if($TreeHash->{$node}{'node_id'} ~~ undef);
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


sub TreeIntersection($$$$){
	
	my ($treeAstring,$treeBstring,$verboseintersection,$TreeAOnly) = @_;
	 #Using IO::String to create io hadles for the newick strings. Do this externally to this function using my $io = IO::String->new($string);
	#$verboseintersection is a flag for printing out a whole load of info regarding the trees and which nodes intersect.
	
	$verboseintersection = 1 if ($verboseintersection ~~ undef);
	
	my ($Aroot,$TreeAObject) = BuildTreeCacheHash($treeAstring);
	my ($Broot,$TreeBObject) = BuildTreeCacheHash($treeBstring);
	
		#Dictionary of node_name to node key
	my $TreeADictionary = {};
	my $TreeBDictionary = {};
	
	
	
	map{$TreeADictionary->{$TreeAObject->{$_}{'node_id'}} = $_;}@{$TreeAObject->{$Aroot}{'Clade_Leaves'}};
	map{$TreeBDictionary->{$TreeBObject->{$_}{'node_id'}} = $_;}@{$TreeBObject->{$Broot}{'Clade_Leaves'}};
	
	my @Ataxa = keys(%$TreeADictionary);
	my @Btaxa = keys(%$TreeBDictionary);
	
	if($verboseintersection){
	
		print  "Tree A leaves: [".scalar(@Ataxa)."] \n";
		print  join(',',sort(@Ataxa));
		print  "\n";
		
		print  "Tree B leaves: [".scalar(@Btaxa)."] \n";
		print  join(',',sort(@Btaxa));
		print  "\n";
	}
	
	my ($Union,$Intersection,$ListAExclusive,$ListBExclusive) = IntUnDiff(\@Ataxa,\@Btaxa);
	
	die "No taxa in commonbetween the two trees!\n" if (scalar(@$Intersection) == 0);

	#A Remove
	foreach my $ANodeToRemove (@$ListAExclusive){

		my $NodeID2Splice = $TreeADictionary->{$ANodeToRemove};
		Splice_Node($TreeAObject,$NodeID2Splice,$Aroot);
	}
	
	sanitise_TreeHash($TreeAObject,$Aroot,$verboseintersection);

unless($TreeAOnly){
	#Bremove
	foreach my $BNodeToRemove (@$ListBExclusive){
		
		my $NodeID2Splice = $TreeBDictionary->{$BNodeToRemove};
		Splice_Node($TreeAObject,$NodeID2Splice,$Aroot);
	}
	
	sanitise_TreeHash($TreeBObject,$Broot,$verboseintersection);
}
	#Sanitise trees
		
	# Output tree descriptions
	
	if($verboseintersection){ #Some useful info regaridng leaves of different trees
		print STDERR "Taxa in common: (".scalar(@$Intersection).")\n";
		my $OutString = join(',',sort(@$Intersection));
		print STDERR $OutString."\n";
		
		print STDERR "Tree A Tree Unique leaves: (".scalar(@$ListAExclusive).")\n";
		$OutString = join(',',@$ListAExclusive);
		print STDERR $OutString."\n";
		
		print STDERR "Output A Tree leaves: (".scalar(@{$TreeAObject->{$Aroot}{'Clade_Leaves'}}).")\n";
		$OutString = join(',',sort(map{$TreeAObject->{$_}{'node_id'}}@{$TreeAObject->{$Aroot}{'Clade_Leaves'}}));
		print STDERR $OutString."\n";
		
		print STDERR "Tree B Tree Unique leaves: (".scalar(@$ListBExclusive).")\n";
		$OutString = join(',',@$ListBExclusive);
		print STDERR $OutString."\n";
		
		print STDERR "Output B Tree leaves: (".scalar(@{$TreeBObject->{$Broot}{'Clade_Leaves'}}).")\n";
		$OutString = join(',',sort(map{$TreeBObject->{$_}{'node_id'}}@{$TreeBObject->{$Broot}{'Clade_Leaves'}}));
		print STDERR $OutString."\n";
	}
	
	return($TreeAObject,$Aroot,$TreeBObject,$Broot);	
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


sub assignLeftRightIDs2TreeHash($$){
	
	my ($TreeCacheHash,$root) = @_;
	
	#Test for binary tree
	
	my $is_binary = isrootedbinary_TreeHash($TreeCacheHash,$root);
	
	if($is_binary){
	
		recursivelyassigneLeftRightids($TreeCacheHash,$root,0);
	
	}else{
	
		die "Tree is non-binary, so you cannot assign left and right ids to it!\n";
	}
	
	return(1); #A list of leaves with the trait of interest, provided in the arguments
}

=pod
=item *assignLeftRightIDs2TreeHash($TreeHash,$root)

Given a treehash, this function will test that it is a rooted binary tree and then assign left and right ids to the tree. 

TODO (Consider a function for non-binary trees at a later date based on natural number indexing)

Creates or reassigns two $TreeHash entries:

$TreeCacheHash->{NodeID}{'left_id'} = node left id
$TreeCacheHash->{NodeID}{'right_id'} = node right id

=cut

sub recursivelyassigneLeftRightids{
	
	my ($TreeCacheHash,$Node,$id_count) = @_;
		
	my @Descendets = @{$TreeCacheHash->{$Node}{'each_Descendent'}};
	
	#Assign left id
	$TreeCacheHash->{$Node}{'left_id'} = $id_count++;
	
	foreach my $Descendent (@Descendets){
		
		$id_count = recursivelyassigneLeftRightids($TreeCacheHash,$Descendent,$id_count);
	}
	
	#Assign_right_id
	$TreeCacheHash->{$Node}{'right_id'} = $id_count++;
	
	return($id_count);
}


=pod
=item *recursivelyassigneLeftRightids($TreeHash,$root)

Don't call directly. This is a recursive function for numbering left and right ids to a binary tree as a treehash. Use assignLeftRightIDs2TreeHash.

=cut

sub isrootedbinary_TreeHash($$){
	
	my ($TreeCacheHash,$root) = @_;
		
	my @NonLeafNodes = grep{! $TreeCacheHash->{$_}{'is_Leaf'}}(@{$TreeCacheHash->{$root}{'each_Descendent'}},$root);
	
	my $is_binary = 1;	
	
	foreach my $InternalNode (@NonLeafNodes){
		
		$is_binary = 0 if(scalar(@{$TreeCacheHash->{$InternalNode}{'each_Descendent'}}) != 2);
	}

	return($is_binary);
}

=pod
=item *isrootedbinary_TreeHash($TreeHash,$root)

Simply tests if a tree stored within a treehash is a rooted binary tree (all nodes have two descendents except for leaves)

=cut


1;

