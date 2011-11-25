#! /usr/bin/perl -w

=head1 NAME

I<.pl>

=head1 USAGE

 .pl [options -v,-d,-h] <ARGS>

=head1 SYNOPSIS

A script to...

=head1 AUTHOR

B<Joe Bloggs> - I<Joe.Bloggs@bristol.ac.uk>

=head1 COPYRIGHT

Copyright 2010 Gough Group, University of Bristol.

=head1 EDIT HISTORY

=cut

# Strict Pragmas
#----------------------------------------------------------------------------------------------------------------
use strict;
use warnings;
#use diagnostics;

# Add Local Library to LibPath
#----------------------------------------------------------------------------------------------------------------
use lib "/home/sardar/bin/perl-libs-custom";


# CPAN Includes
#----------------------------------------------------------------------------------------------------------------
=head1 DEPENDANCY
B<Getopt::Long> Used to parse command line options.
B<Pod::Usage> Used for usage and help output.
B<Data::Dumper> Used for debug output.
=cut
use Getopt::Long;                     #Deal with command line options
use Pod::Usage;                       #Print a usage man page from the POD comments after __END__
use Data::Dumper;                     #Allow easy print dumps of datastructures for debugging
#use XML::Simple qw(:strict);          #Load a config file from the local directory
use DBI;
use Supfam::Utils;

use TraP::Cluster;

# Command Line Options
#----------------------------------------------------------------------------------------------------------------

my $verbose; #Flag for verbose output from command line opts
my $debug;   #As above for debug
my $help;    #Same again but this time should we output the POD man page defined after __END__
my $file ;
my $method;

#Set command line flags and parameters.
GetOptions("verbose|v!"  => \$verbose,
           "debug|d!"  => \$debug,
           "file|f=s"  => \$file,
           "help|h!" => \$help,
            "method|m=s" => \$method,
        ) or die "Fatal Error: Problem parsing command-line ".$!;



#Print out some help if it was asked for or if no arguments were given.
#pod2usage(-exitstatus => 0, -verbose => 2) if not $file or $help;

# Sub definitions
#----------------------------------------------------------------------------------------------------------------
=head1 DESCRIPTION

Detailed info about the script goes here

=head2 Methods
=over 4
=cut

=item * func
Function to do something
=cut
sub func {
	return 1;
}

# Main Script Content
#----------------------------------------------------------------------------------------------------------------


##Using the viral dataset - lets see how it does

open DATAFILE,"<$file" or die $!;

my $firstline = <DATAFILE>;

my $RawDataHash={};

my @DataLabels = split(/\s/,$firstline);

$firstline = <DATAFILE>; #Skip through another pointless line

while (my $line = <DATAFILE>){
	
	my @SplitLine = split(/\t+/,$line);
	
	my $SpeciesName = shift(@SplitLine);
	
	$RawDataHash->{$SpeciesName}=\@SplitLine;
}

my ($ClusterPositionsHash,$XYClusterGroups) = SOMcluster($RawDataHash,$method,0);

EasyDump("./ClusterPositionsHash.dat", $ClusterPositionsHash);
EasyDump("./XYClusterGroups.dat", $XYClusterGroups);

__END__

