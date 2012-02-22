package Supfam::Utils;
require Exporter;

=head1 NAME

Supfam::Utils

=head1 SYNOPSIS

A few useful functions for scripting, like easier dumping and reading in of data
use Supfam::Utils;

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
					EasyDump
					EasyUnDump
					IntUnDiff
					TabSepFile
					choose_weighted
					CommaSepFile
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use lib "$ENV{HOME}/bin/perl-libs-custom";

use Data::Dumper;
use Term::ProgressBar;
use Math::Combinatorics;
use Supfam::SQLFunc;
use Carp qw(croak);
use Params::Validate qw(:all);


sub EasyDump($$){
	my ($FileName,$Pointer) = @_;
	
	open FH, ">".$FileName or die "Unable to initalise File Handle";
	print FH Dumper($Pointer);
	close FH;
}

=pod
=item * EasyDump(FileName, Pointer)
A wrapper function for the dumper module. Pass in the desired file name and a pointer to the object to be outputted. Does not return.
=cut

sub EasyUnDump($){
	my ($FileName) = shift();
	
	open FH , "<".$FileName or die "Couldn't open $FileName: $!";
	my $FileDump;
	
		while (my $line = <FH>){
			$FileDump .= $line;
		}
	close FH;
	
	my $VAR1;
	eval($FileDump);
	return($VAR1)
}

=pod
=item * EasyUnDump($)
A wrapper function for the dumper module. Given a file produced using dumper, this function will return a pointer to the structure in memory.
=cut

sub IntUnDiff($$){
	
	my ($ListA,$ListB) = @_;
	
	my $UnionHash = {};
	my $Union = [];
	my $Intersection = [];
	my $ListAExclusive = [];
	my $ListBExclusive = [];
	
	my $ListALookup={};
	my $ListBLookup={};
	
	map{$ListALookup->{$_} = 1}@$ListA;
	map{$ListBLookup->{$_} = 1}@$ListB;
	#Initialise a hash for lookup later
	
	foreach my $element (keys(%$ListALookup), keys(%$ListBLookup)) { $UnionHash->{$element}++; }
	#USe the keys of the above hash to get the unique elements in that set
	
	@$Union = keys(%$UnionHash);
	
	 foreach my $element (@$Union) {
	 	   
	 	 if ($UnionHash->{$element} == 2) {  #i.e if it's in both sets      
	 	 	 push (@$Intersection, $element);
	 	 } else {
	 	 
	 	 	if (exists($ListALookup->{$element})){
	 	 		push(@$ListAExclusive, $element);   
	 	 	}elsif(exists($ListBLookup->{$element})){
	 	 		push(@$ListBExclusive, $element); 
	 	 	}else{
	 	 		die "Error in code!\n";
	 	 	}
	 	 }
	 }
	
	return($Union,$Intersection,$ListAExclusive,$ListBExclusive);

}

=pod
=item * IntUnDiff($$)
A quick function to calculate some basic set statistics between two lists (supplied as pointers to two arrays in). Returns four
pointers to arrays of the 1. Union 2. Intersection 3. Elements unique to list A 4. Elements unique to list B.

NOTE - repeats in input list are ignored! So the output lists are unique! so if list A = (A,A,B,C) and list B = (A,D,E,F), union = (A,B,C,D,E,F) and intersection (A) etc. 

=cut


sub TabSepFile{
	
	my ($Fields,$OutputData,$fileout,$totals,$defualtvalue) = @_;

	$defualtvalue = 'N/A' unless defined($defualtvalue);
	
	# $fields = [field headings]
	# $Output data = {row titles => {field1 => val1, field2 => va2 ...}}. For fields not included in the OutputData hash, $defualtvalue will be assumed to be the value
	#$totals = {field1 = total_val1, field2 => total_val2 ...}
	
	open TABOUTPUT, ">$fileout";
	
	print TABOUTPUT "Entry\t";
	print TABOUTPUT join("\t",@$Fields);
	print TABOUTPUT "\n";
	
	if(defined($totals)){
		
		my @FieldValues = @{$totals}{@$Fields}; #Only retrieve the values for which there are fields in @$Fields
		
		print TABOUTPUT "Field_Sum\t";
		print TABOUTPUT join("\t",@FieldValues);
		print TABOUTPUT "\n";
	}
	
	my $FieldsHash = {};
	map{$FieldsHash->{$_}=$defualtvalue}@$Fields;
	#Initialise a hash with default values
	
	foreach my $Entry (keys(%$OutputData)){
	
		my $EntrySpecificFieldHash ={};
		%$EntrySpecificFieldHash = %$FieldsHash ;
		
		@{$EntrySpecificFieldHash}{keys(%{$OutputData->{$Entry}})}=values(%{$OutputData->{$Entry}}); #update entry specific hash using a hash slice
		
		my @FieldValues = @{$EntrySpecificFieldHash}{@$Fields}; #Only retrieve the values for which there are fields in @$Fields

		print TABOUTPUT "$Entry\t";
		print TABOUTPUT join("\t",@FieldValues);
		print TABOUTPUT "\n";
	}
	
	close TABOUTPUT;

	return(1);
	
}

=pod
=item * TabSepFile($Fields,$OutputData,$fileout,$totals -optional,$defualtvalue - optional)

An easy way to dump a whole load of Entry=>{field1 =>val1, field2 => val2 ...} data to a tab seperated file. First line of output is a list of fields (as specieifed in @$Fields), followed by a line per key in
%$OutputData. Only fields in @$Fields are outputted. $Output data = {row titles => {field1 => val1, field2 => va2 ...}}. For fields not included in the OutputData hash, $defualtvalue will be assumed to be the value

File save to $fileout

=cut

sub CommaSepFile{
	
	my ($Fields,$OutputData,$fileout,$totals,$defualtvalue) = @_;
	
	$defualtvalue = 'N/A' unless defined($defualtvalue);
	
	# $fields = [field headings]
	# $Output data = {row titles => {field1 => val1, field2 => va2 ...}}. For fields not included in the OutputData hash, $defualtvalue will be assumed to be the value
	#$totals = {field1 = total_val1, field2 => total_val2 ...}
	
	
	open COMMAOUTPUT, ">$fileout";
		
	print COMMAOUTPUT "Entry,";
	print COMMAOUTPUT join(",",@$Fields);
	print COMMAOUTPUT "\n";
	
	
	if(defined($totals)){
		
		my @FieldValues = @{$totals}{@$Fields}; #Only retrieve the values for which there are fields in @$Fields
		
		print COMMAOUTPUT "Field_Sum,";
		print COMMAOUTPUT join(",",@FieldValues);
		print COMMAOUTPUT "\n";
	}
	
	
	my $FieldsHash = {};
	map{$FieldsHash->{$_}=$defualtvalue}@$Fields;
	
	foreach my $Entry (keys(%$OutputData)){
	
		my $EntrySpecificFieldHash ={};
		%$EntrySpecificFieldHash = %$FieldsHash ;
		
		@{$EntrySpecificFieldHash}{keys(%{$OutputData->{$Entry}})}=values(%{$OutputData->{$Entry}}); #update entry specific hash using a hash slice
		
		my @FieldValues = @{$EntrySpecificFieldHash}{@$Fields}; #Only retrieve the values for which there are fields in @$Fields
		
		print COMMAOUTPUT "$Entry,";
		print COMMAOUTPUT join(",",@FieldValues);
		print COMMAOUTPUT "\n";
	}
	
	close COMMAOUTPUT;
}

=pod
=item * CommaSepFile($Fields,$OutputData,$fileout,$totals -optional,$defualtvalue - optional)

An easy way to dump a whole load of Entry=>{field1 =>val1, field2 => val2 ...} data to a c seperated file. First line of output is a list of fields (as specieifed in @$Fields), followed by a line per key in
%$OutputData. Only fields in @$Fields are outputted. $Output data = {row titles => {field1 => val1, field2 => va2 ...}}. For fields not included in the OutputData hash, $defualtvalue will be assumed to be the value

File save to $fileout

=cut


sub choose_weighted{
     validate_pos(@_, 
		  { type => ARRAYREF },
		  { type => CODEREF | ARRAYREF}
 	);

    my ($objects, $weightsArg ) = @_;
    my $calcWeight = $weightsArg if 'CODE' eq ref $weightsArg;
    my @weights;		# fix wasteful of memory
    if( $calcWeight){
	@weights =  map { $calcWeight->($_) } @$objects; 
    }
    else{
	@weights =@$weightsArg;
	if ( @$objects != @weights ){
	    croak "given arefs of unequal lengths!";
	}
    }

    my @ranges = ();		# actually upper bounds on ranges
    my $left = 0;
    for my $weight( @weights){
	$weight = 0 if $weight < 0; # the world is hostile...
	my $right = $left+$weight;
	push @ranges, $right;
	$left = $right;
    }
    my $weightIndex = rand $left;
    for( my $i =0; $i< @$objects; $i++){
	my $range = $ranges[$i];
	return $objects->[$i] if $weightIndex < $range;
    }
}

=pod
=item * choose_weighted ($object_Aref, $weights_Aref )

Just one function, a simple means of making a weighted random choice

The implementation uses rand to calculate random numbers.

	choose_weighted ($object_Aref, $weights_Aref )
or 
	choose_weighted ($object_Aref, $weights_codeRef )

In the second case, the coderef is called on each object to determine its weight;

=cut


1;

