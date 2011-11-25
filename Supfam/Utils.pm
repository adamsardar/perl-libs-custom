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
					CommaSepFile
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use Data::Dumper;
use Term::ProgressBar;
use Math::Combinatorics;
use Supfam::SQLFunc;


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
	
	my $switch = 0;
	#A flag for if the two lists are flipped around in the next step. This allows for correct reporting of the return values
	
	if (scalar(@$ListA) > scalar(@$ListB)){
		my $TempList = $ListA;
		$ListA = $ListB;
		$ListB = $TempList;
		$switch=1;
	}
	# This is to make sure that the code runs efficiently, which requires ListA to be the smaller of the two lists.
	
	my $UnionHash = {};
	my $Union = [];
	my $Intersection = [];
	my $ListAExclusive = [];
	my $ListBExclusive = [];
	
	my $ListALookup={};
	
	foreach(@$ListA){$ListALookup->{$_} = 1;}
	#Initialise a hash for lookup later
	
	foreach my $element (@$ListA, @$ListB) { $UnionHash->{$element}++; } 
	
	@$Union = keys(%$UnionHash);
	
	 foreach my $element (@$Union) {     
	 	   
	 	 if ($UnionHash->{$element} == 2) {  #i.e if it's in both sets      
	 	 	 push (@$Intersection, $element);     
	 	 } else {     
	 	 	no warnings 'uninitialized';
	 	 	#This is to stop Perl moaning about elements not beining initialised in the lookup hash below
	 	 
	 	 	if ($ListALookup->{$element}){
	 	 		push(@$ListAExclusive, $element);   
	 	 	}else {
	 	 		push(@$ListBExclusive, $element); 
	 	 	}
	 	 } 
	 }
	
	unless ($switch) {return($Union,$Intersection,$ListAExclusive,$ListBExclusive);
	}else {			return($Union,$Intersection,$ListBExclusive,$ListAExclusive);
	};
}

=pod
=item * IntUnDiff($$)
A quick function to calculate some basic set statistics between two lists (supplied as pointers to two arrays in). Returns four
pointers to arrays of the 1. Union 2. Intersection 3. Elements unique to list A 4. Elements uniqur to list B.
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

1;

