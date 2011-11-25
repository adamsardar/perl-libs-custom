package Supfam::SQLFunc;
require Exporter;

=head1 NAME

Supfam::SQLFunc.pm

=head1 SYNOPSIS

Holds all the functions required to do interesting things with the Superfamily database.
use Supfam::SQLFunc;

=head1 AUTHOR

Matt Oates (Matt.Oates@bristol.ac.uk)

=head1 COPYRIGHT

Copyright 2010 Gough Group, University of Bristol.

=head1 SEE ALSO

Supfam::Config.pm

=head1 DESCRIPTION

=cut

our @ISA       = qw(Exporter AutoLoader);
our @EXPORT    = qw(
			dbConnect
			dbDisconnect
                  );
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;

use DBI;
use Data::Dumper;
use Term::ProgressBar;
use Math::Combinatorics;

use Supfam::Config;

=pod
=head2 Methods
=over 4
=cut


sub dbConnect {
	
return DBI->connect("DBI:mysql:dbname=$SUPFAM{'database'}{'name'};host=$SUPFAM{'database'}{'host'}"
                                        ,$SUPFAM{'database'}{'user'}
                                        , undef
                                        ,{RaiseError =>1}
                                    ) or die "Fatal Error: couldn't connect to $SUPFAM{'database'}{'name'} on $SUPFAM{'database'}{'host'}";
}


=pod

=back

=cut

sub dbDisconnect {
	my $dbh = shift;
	#warn "Closing Database Connection!\n";
	return $dbh->disconnect();
}



1;

