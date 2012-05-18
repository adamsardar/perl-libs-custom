#!/usr/bin/env perl

package Package::Skeleton;

require Exporter;

our @ISA = qw(Exporter);

our %EXPORT_TAGS = (
'all' => [ qw(
			sub1
			sub2
) ],
'yourtag' => [ qw(
			sub1
) ],
);
our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
our @EXPORT = qw();

our $VERSION   = 1.00;

use strict;
use warnings;

=head1 NAME

Pacakge::Skeleton v1.0 - Skeleton module for making new perl modules is sardar repositories

=head1 DESCRIPTION

Just a skeleton layout for each module to start from.

=head1 EXAMPLES

use Pacakge::Skeleton qw/all/;

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



=head1 FUNCTIONS DEFINED

=over 4
=cut

=item * sub1
Function to do something
=cut
sub sub1 {
    my ($var) = @_;
	return 1;
}

=item * sub2
Function to do something
=cut
sub sub2 {
    my ($var) = @_;
	return 2;
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

