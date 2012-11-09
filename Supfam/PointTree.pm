#!/usr/bin/env perl

package Supfam::PointTree;

our $VERSION   = 1.00;

use strict;
use warnings;

=head1 NAME

Supfam::PointTree

=head1 DESCRIPTION

An implementation of point trees in perl. The core reason to do this is to replace looping through indicies of an array (which is O(n)) with
 a balanced tree search (order O(log(n))). Note that tree construction takes O(nlog{n}), so this is oly useful in highly repetative loop regions.
 
=head1 EXAMPLES

use Supfam::PointTree qw/all/;

=head1 AUTHOR

B<Adam Sardar> - I<adam@sardar.me.uk>

=head1 NOTICE

B<Adam Sardar> (2012) First features added.

=head1 LICENSE AND COPYRIGHT

B<Copyright 2012 Adam Sardar>

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

#use lib '~/lib';

=head1 DEPENDANCY

B<Data::Dumper> Used for debug output.

=cut
use Data::Dumper; #Allow easy print dumps of datastructures for debugging
use Carp;
use List::Util;
use POSIX ('ceil');
use Math::Random ('random_uniform');

our $AUTOLOAD;  # it's a package global


my %fields = (
    TREE        => undef,
    LOCK         => 0,
    NPOINTS       => 0,
    DELPOINTS => [],
    NUMBEROFDELPOINTS => 0,
);


=head1 FUNCTIONS DEFINED

=over 4
=cut

=item * new
The constructor method for the class Supfam::PointTree
=cut
sub new {
		
		my $proto = shift;
        my $class = ref($proto) || $proto;
        
        my $self  = {
        	 "_permitted" => \%fields,
            %fields,
        };

        bless ($self, $class);
        return($self);
}

=item * AUTOLOAD
AUTOLOADER for Supfam::PointTree
=cut

sub AUTOLOAD {
	
	my $self = shift;
	my $type = ref($self) or croak "$self is not an object";
	
	(my $name = $AUTOLOAD) =~ s/.*://;   # strip fully-qualified portion
	unless (exists $self->{"_permitted"}{$name} ) {
		croak "Can't access `$name' field in class $type";	
	}
	
	if (@_) {
		return $self->{$name} = shift;
	} else {
		return $self->{$name};
	}
}



=item * build
Method to add intervals into our point tree
=cut

sub build {
	
	my $self = shift;
	my $Intervals = shift;
	
	croak "You may only build a point tree once. Addition of intervals is not supported at this stage - add them all in one go.\n" if($self->LOCK);
	
	#Build the tree and place it into $self->TREE
	
	my $Tree = {};
	enter_tree($Tree,$Intervals);
	#Physically build a blanced binary tree of intervals.
	
	$self->TREE($Tree);
	$self->NPOINTS(scalar(@$Intervals));
	$self->LOCK(1);
	$self->DELPOINTS([]);
	
	return 1;
}

=item * enter_tree
A sub-routine called by the build method so as to actually construct a tree representation of the intervals. Please do not call call directl (it is not exported by this package and is not a method).
Works through a depth-first recursion.
=cut

sub enter_tree{
	
	carp "This sub takes 2 arguments exactly\n" unless(scalar(@_) == 2);

	my ($Tree,$PointsOnLevel) = @_;
	# Tree is a hash form of the tree under construction, $PointsOnLevel is an array ref of [list of items on this level]
	
	unless(scalar(@$PointsOnLevel) > 0){
		
		carp "No points on this level!\n";
	}
	
	
	$Tree->{'PointsOnLevel'}= scalar(@$PointsOnLevel);
	
	if($Tree->{'PointsOnLevel'} == 1){
				
			$Tree->{'Item'} = $$PointsOnLevel[0];
			
			return(1);
		
	}elsif($Tree->{'PointsOnLevel'} > 1){
				
		#Sort items
		my @SortedPoints = sort{$a <=> $b}@$PointsOnLevel;
		
		#Split points into left and right lineages
		my ($LeftLineage,$RightLineage) = ([],[]);
					
		#Set Divinding Point between two groups
		my $Midpoint = ceil($Tree->{'PointsOnLevel'}/2);
		my $MidpointIndex = $Midpoint-1;
		
		$Tree->{'Dividing_Point'} = $SortedPoints[$MidpointIndex];

		#Collect Left Points of divider
		@$LeftLineage=@SortedPoints[0 .. $MidpointIndex];
		
		#Right Points of divider
		@$RightLineage=@SortedPoints[$MidpointIndex+1 .. scalar(@SortedPoints)-1];
		
		my ($LeftSubTree,$RightSubTree) = ({},{});
		
		$Tree->{'LeftTree'} = $LeftSubTree;
		$Tree->{'RightTree'} = $RightSubTree;
		
		#Left Tree Lineage (note, recursive function)
		enter_tree($LeftSubTree,$LeftLineage);
		
		#Right Tree Lineage (note, recursive function)
		enter_tree($RightSubTree,$RightLineage);

	}else{
		
		croak "Somethign seriously wrong with point tree construction - less than one interval on this level!\n";
	}
	
}


=item * Search
Given an array ref of points on an interval used to construct the PointTree, this function will find them all and return the values back in the arrayref provided ($IntervalPoints)
=cut
sub Search {
   
    my $self = shift;
    my $PointsToSearch = shift;
    my $IntervalPoints = shift;
    
    carp "This method takes a single pointer to an array as function input. Recieved a ".ref($PointsToSearch)."\n" unless (ref($PointsToSearch) eq 'ARRAY');
    
    #Search to look through tree and find point of interest
   
    my $Tree = $self->TREE;

	while(my $PointToSearch = pop(@$PointsToSearch)){

		my $IntervalPoint; #This will be the returned value. We propogate a pointer to this down through the levels of the following tree search
	
		Point_Search_Iteratively($Tree,$PointToSearch,\$IntervalPoint);
	    
	    #print $IntervalPoint;
	 	push(@$IntervalPoints,$IntervalPoint);
	}
}

=item * Point_Search_Recursively
A sub called by the Search method. Please don't call manually. This will seek the interval/divider immediitaly above the point given - this is useful for mapping, say, continuous points
along a given interval to discrete times (e.g. deletion events to nodes on a tree).
=cut

sub Point_Search_Recursively{
   
   carp "This sub takes 3 arguments exactly\n" unless(scalar(@_) == 3);
   
    my ($Tree,$Point,$IntervalPointPointer) = @_;
	#We propagate $IntervalPointPointer down through the depth first search of the point tree. We then overwrite that memory location and then collapse back up the tree.
	
	if($Tree->{'PointsOnLevel'} == 1){
		
		$$IntervalPointPointer = $Tree->{'Item'}; 
		#Collect the point on this level and write to the memory location of $IntervalPointPointer
		
	}else{
		
		my $SubTreeToSearch;
		
		if($Point > $Tree->{'Dividing_Point'}){
			
			#Right Side
			$SubTreeToSearch =$Tree->{'RightTree'};
			
		}else{
			
			#Left Side
			$SubTreeToSearch =$Tree->{'LeftTree'};
		}
		
		Point_Search_Recursively($SubTreeToSearch,$Point,$IntervalPointPointer);
	}
}


=item * Point_Search_Iteratively
A sub called by the Search method. Please don't call manually. This will seek the interval/divider immediitaly above the point given - this is useful for mapping, say, continuous points
along a given interval to discrete times (e.g. deletion events to nodes on a tree).
=cut


sub UniformAssign{
	
	my $self = shift;
    my $NumberOfPoints = shift;
    
    my $UniformDeletions = [];
    @$UniformDeletions = random_uniform($NumberOfPoints,0,1);
 
	my $DeletionPoints = $self->DELPOINTS;
	$self->NUMBEROFDELPOINTS($NumberOfPoints);

	$self->Search($UniformDeletions,$DeletionPoints);#Find the node directly below the deletion
}

=item * UniformAssign
Populates a stock array of deletion points across a subtree that can be drawn from using UniformDraw. This is done so that a pool of deletions can be constructed, and then simply pulled from as required.
=cut

sub UniformDraw{
	
	my $self = shift;
    my $NumberOfDraws = shift;
	my $ReturnedDelPointsArray = shift;

    
    my $DelPointsArray = $self->{'DELPOINTS'};
    #Just to avoid AUTOLOAD and its speed drop. This is bad, bad, bad coding form. But perl is shite for OO, so what are you gonna do?
    $self->UniformAssign($self->NUMBEROFDELPOINTS) if(scalar(@$DelPointsArray) < $NumberOfDraws);

	my $count = 0;
	
	while ($count < $NumberOfDraws){
		
		push(@$ReturnedDelPointsArray,pop($DelPointsArray));	
		$count++;
	}

}

=item * UniformDraw
Pull a whole load of uniform dleetions from a pre-stored pool, then map them back to the tree.
=cut

sub Point_Search_Iteratively{
   
   carp "This sub takes 3 arguments exactly\n" unless(scalar(@_) == 3);
   
    my ($Tree,$Point,$IntervalPointPointer) = @_;
	#We propagate $IntervalPointPointer down through the depth first search of the point tree. We then overwrite that memory location and then collapse back up the tree.
	
	while($Tree->{'PointsOnLevel'} != 1){
		
		if($Point > $Tree->{'Dividing_Point'}){
			
			#Right Side
			$Tree =$Tree->{'RightTree'};
			
		}else{
			
			#Left Side
			$Tree =$Tree->{'LeftTree'};
		}
	}
	
	$$IntervalPointPointer = $Tree->{'Item'}; 
}


=item * DESTROY
The destructor method for the class Supfam::IntervalTree
=cut

sub DESTROY {

}

=pod

=back

=head1 TODO

=over 4

=item Add feature here...

=back

=cut

1;
__END__

