package TraP::Cluster;
require Exporter;

=head1 NAME

TraP::Cluster

=head1 SYNOPSIS

First perl module in a library that will become that underpinning the clustering portion of the TraP database

use TraP::Cluster;

=head1 AUTHOR

Adam Sardar (adam.sardar@bris.ac.uk)
Owen Rackham (owen.rackham@bris.ac.uk)
Matt Oates (matt.oates@bris.ac.uk)
David Morais (david.morais@bris.ac.uk)

=head1 COPYRIGHT

Copyright 2011 Gough Group, University of Bristol.

=head1 SEE ALSO


=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
					SOMcluster
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use lib "/home/sardar/bin/perl-libs-custom/";

use POSIX;

use Supfam::SQLFunc;
use Algorithm::Cluster;
use Supfam::Utils;

sub SOMcluster($$$){
	
	my ($RawDataHash,$ClusterDistMethod,$transpose) = @_;
	
	# $RawdataHash->{key}=[values],  $ClusterDistMethod is the distance method to use in the self organising map, $Transpose is a flag for whether to invert the data matrix
	
	my @AllowedClusterMethods = qw(c a u x s k e b);
	#c: correlation a: absolute value of the correlation u: uncentered correlation x: absolute uncentered correlation s: Spearman''s rank correlation k: Kendall''s tau e: Euclidean distance b: City-block distance
	die "Incorrect clustering method chosen\n" unless (grep{m/$ClusterDistMethod/}@AllowedClusterMethods);
	
	my @DataLabels = keys(%$RawDataHash);
	my $NoOfdataPoints = scalar(@DataLabels);
		
	my $nxgrid = ceil(sqrt($NoOfdataPoints)); #http://www.viscovery.net/faqs suggests that the number of nodes in a SOM should be equal to the number of datapoints divided by ten
	my $nygrid = $nxgrid; #Creating a square initial SOM grid
	
	#Create 2D data array with labels
	
	my $DataValues2Cluster = [(@{$RawDataHash}{@DataLabels})]; #Creates a pointer to a 2D array of data values, in the order specieifed by @DataLabels
	
	my %param = (
				  'data' =>  $DataValues2Cluster,
				  'transpose' => $transpose, 
				  'nxgrid' => $nxgrid,
				  'nygrid' => $nygrid,
				  'niter' => 100,
				  'dist' => "$ClusterDistMethod"
	);
	
	my ($clusterid) = Algorithm::Cluster::somcluster(%param);
	#clusterid is of form [[x1,y1],[x2,y2],...,[xn,yn]]
	
	my $ClusterPositionsHash = {};
	@{$ClusterPositionsHash}{@DataLabels}=@$clusterid;
	#Output in form $ClusterPositionsHash->{data_label}=[x,y] where x and y are the coordinates of the data point on the SOM

	my $ReducedClusterPositionsHash = {};
	map{$ReducedClusterPositionsHash->{$_}=join(',',@{$ClusterPositionsHash->{$_}})}@DataLabels;
	#Very similar to the above $ClusterPositionsHash, but now the value is not an array ref but a string "x,y"
	
	my $XYClusterGroups ={};
	
	foreach my $DataLabel (@DataLabels){
		
		my $XYPoint = $ReducedClusterPositionsHash->{$DataLabel};
		
		$XYClusterGroups->{$XYPoint}=[] unless(exists($XYClusterGroups->{$XYPoint}));
		push(@{$XYClusterGroups->{$XYPoint}},$DataLabel);
	}
	#Taking the above $ClusterPositionsHash, $XYClusterGroups will be of the form $XYClusterGroups->{"x,y"}=[samples] - i.e. it inverts $ReducedClusterPositionsHash
	
	return($ClusterPositionsHash,$XYClusterGroups);
}

=pod
=item * SOMcluster($RawDataHash,$ClusterDistMethod,$transpose)

Input:
$RawdataHash->{key}=[values],  $ClusterDistMethod is the distance method to use in the self organising map, $Transpose is a flag for whether to invert the data matrix

This is really just a wrapper to the Algorithm::Cluster perl module but with data input and output in the form wanted for TraP projects.

Returns: $ClusterPositionsHash,$XYClusterGroups

$ClsuterPositionsHash->{data_label}=[x,y]
$XYClusterGroups->{"x,y"}=[data_labels] i.e. a stringified x,y coord to data lables 
=cut

1;

