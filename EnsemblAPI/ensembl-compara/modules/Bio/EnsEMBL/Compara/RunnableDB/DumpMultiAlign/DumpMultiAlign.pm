=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::DumpMultiAlign::DumpMultiAlign

=head1 SYNOPSIS

This RunnableDB module is part of the DumpMultiAlign pipeline.

=head1 DESCRIPTION

This RunnableDB module runs DumpMultiAlign jobs. It creates emf2maf jobs if
necessary and compression jobs

=cut

package Bio::EnsEMBL::Compara::RunnableDB::DumpMultiAlign::DumpMultiAlign;

use strict;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;
use base ('Bio::EnsEMBL::Compara::RunnableDB::BaseRunnable');

=head2 strict_hash_format

    Description : Implements strict_hash_format() interface method of Bio::EnsEMBL::Hive::Process that is used to set the strictness level of the parameters' parser.
                  Here we return 0 in order to indicate that neither input_id() nor parameters() is required to contain a hash.

=cut

sub strict_hash_format {
    return 0;
}

sub fetch_input {
    my $self = shift;
}

sub run {
    my $self = shift;

    my $cmd = $self->param('cmd');

    #append full path to output_file
    my $full_output_file = $self->param('output_dir') . "/" . $self->param('output_file');
    $cmd .= " --output_file $full_output_file";

    #Write a temporary file to store gabs to dump
    if ($self->param('start') && $self->param('end')) {
	$self->_write_gab_file();
	$cmd .= " --file_of_genomic_align_block_ids " . $self->param('tmp_file');
    }

    #substitute any hashes in analysis parameters with the correct values from analysis_job
    $cmd = $self->param_substitute($cmd);
    #print "cmd $cmd \n";

    #
    #Run DumpMultiAlign cmd
    #
    if(my $return_value = system($cmd)) {
        $return_value >>= 8;
        die "system( $cmd ) failed: $return_value";
    }
    #
    #Check number of genomic_align_blocks written is correct
    # 
    $self->_healthcheck();
}

sub write_output {
    my $self = shift @_;

    #delete tmp file
    unlink($self->param('tmp_file'));

    #
    #Create emf2maf job if necesary
    #
    if ($self->param('maf_output_dir')) {
	my $output_ids = "{\"output_file\"=>\"" . $self->param('dumped_output_file') . "\", \"num_blocks\" => \"" . $self->param('num_blocks') . "\"}";

	$self->dataflow_output_id($output_ids, 2);

    } else {
	#Send dummy jobs to emf2maf
	$self->dataflow_output_id("{}", 2);
    }

    #
    #Create Compress jobs
    #
    my $output_ids = "{\"output_file\"=>\"" . $self->param('dumped_output_file') . "\"}";
    $self->dataflow_output_id($output_ids, 1);
}

#
#Check the number of genomic_align_blocks written is correct
#
sub _healthcheck {
    my ($self) = @_;
    
    #Find out if split into several files
    my $dump_cmd = $self->param('extra_args');
    my $chunk_num = $dump_cmd =~ /chunk_num/;
    my $output_file = $self->param('output_dir') . "/" . $self->param('output_file');

    #not split by chunk eg supercontigs so need to check all supercontig* files
    if (!$chunk_num) {
	if ($output_file =~ /\.[^\.]+$/) {
	    $output_file =~ s/(\.[^\.]+)$/_*$1/;
	}
    } else {
	#Have chunk number in filename
	$output_file = $self->param('output_dir') . "/" . $self->param('dumped_output_file');
    }

    my $cmd;
    if ($self->param('format') eq "emf") {
	$cmd = "grep DATA " . $output_file . " | wc -l";

    } elsif ($self->param('format') eq "maf") {
	$cmd = "grep ^a " . $output_file . " | wc -l";
    }
    my $num_blocks = `$cmd`;
    chomp $num_blocks;
    if ($num_blocks != $self->param('num_blocks')) {
	die("Number of block dumped is $num_blocks but should be " . $self->param('num_blocks'));
    } else {
	print "Wrote " . $self->param('num_blocks') . " blocks\n";
	#Store results in table. Not really necessary but good to have 
	#visual confirmation all is well
	my $sql = "INSERT INTO healthcheck (filename, expected,dumped) VALUES (?,?,?)";
	my $sth = $self->analysis->adaptor->dbc->prepare($sql);
	$sth->execute($self->param('output_file'), $self->param('num_blocks'), $num_blocks);
	$sth->finish();
    }
}

#
#Write temporary file containing a list of genomic_align_block_ids for 
#inputting into DumpMultiAlign
#
sub _write_gab_file {
    my ($self) = @_;

    my $sql = "SELECT * FROM other_gab WHERE genomic_align_block_id BETWEEN ? AND ?";
    my $sth = $self->analysis->adaptor->dbc->prepare($sql);
    $sth->execute($self->param('start'), $self->param('end'));
    
    my $tmp_file = "/tmp/other_gab_$$.out";
    $self->param('tmp_file', $tmp_file);
    
    open(FILE, ">$tmp_file") || die ("Couldn't open $tmp_file for writing"); 

    while (my $row = $sth->fetchrow_arrayref) {
	print FILE $row->[0] . "\n";

    }
    close(FILE);
    $sth->finish;
}

1;
