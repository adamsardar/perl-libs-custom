#!/usr/bin/env perl

# This is a program that takes a file with a list of stable IDs (not
# exon stable IDs), and outputs a comma-separated list of the history of
# each of these stable IDs.  The history ends with either the current
# release or at the point when the stable ID was retired.

use strict;
use warnings;

use lib "/home/sardar/bin/perl-libs-custom/EnsemblAPI/";
use lib "/home/sardar/bin/perl-libs-custom/EnsemblAPI/ensembl/modules";

use Bio::EnsEMBL::Registry;
use IO::File;
use Getopt::Long;
use DBI qw(:sql_types);

my ( $filename, $species );
my $help = '';

if ( !GetOptions( 'file|f=s'    => \$filename,
                  'species|s=s' => \$species,
                  'help|h!'     => \$help )
     || !defined($species)
     || $help )
{
  print <<END_USAGE;

Usage:
  $0 --species=species [ --file=filename ] [ --species=species ]

  $0 --help


    --species / -s  Name of species.

    --file    / -f  (Optional) Name of file containing a list of stable
                    IDs (new line seperated).  The default is to read the
		    list from standard  input.

    --help    / -h  To see this text.

Example usage:

  $0 -s mouse -f idlist.txt

  # Same as the above:
  $0 -s mouse < idlist.txt

  sort -u longlist.txt | $0 -s human

  Output is printed to STDOUT (*sigh*) as a comma-separated list of the history of
  each of these stable IDs.  The history ends with either the current
  release or at the point when the stable ID was retired.

END_USAGE

  exit(1);
} ## end if ( !GetOptions( 'file|f=s'...))

$filename ||= '-';

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db( -host => 'ensembldb.ensembl.org',
                                  -user => 'anonymous' );

my $adaptor = $registry->get_DBAdaptor( $species, 'Core' );

my $in = IO::File->new($filename);

if ( !defined($in) ) {
  die( sprintf( "Could not open file '%s' for reading", $filename ) );
}

# We could do what we want to do with the API, but this is simpler and
# quicker, at the moment.  As always, when using plain SQL against our
# databases, the user should not be surprised to see the code break when
# we update the schema...
my $statement = q(
SELECT  old_version,
        old_release,
        new_stable_id, new_version,
        new_release,
        score
FROM    stable_id_event
  JOIN  mapping_session USING (mapping_session_id)
WHERE   old_stable_id = ?
ORDER BY old_version ASC, CAST(new_release AS UNSIGNED)
);

my $sth = $adaptor->dbc()->db_handle()->prepare($statement);

while ( my $stable_id = $in->getline() ) {
  chomp($stable_id);

  # Strip off any comment (from '#' to the end of the line).
  $stable_id =~ s/\s*#.*$//;

  # Skip lines containing only whitespace.
  if ( $stable_id =~ /^\s*$/ ) { next }

  print("Old stable ID, New stable ID, Release, Mapping score\n");

  $sth->bind_param( 1, $stable_id, SQL_VARCHAR );

  $sth->execute();

  my ( $version, $release, $new_stable_id, $new_version, $new_release,
       $score );

  $sth->bind_columns( \( $version,     $release,     $new_stable_id,
                         $new_version, $new_release, $score ) );

  while ( $sth->fetch() ) {
    if ( defined($new_stable_id) ) {
      printf( "%s.%s, %s.%s, %s, %s\n",
              $stable_id, $version, $new_stable_id, $new_version,
              $new_release, $score );
    } elsif ( !defined($new_stable_id) ) {
      printf( "%s.%s, <retired>, %s, %s\n",
              $stable_id, $version, $new_release, $score );
    }
  }

  print("\n");
} ## end while ( my $stable_id = $in...)
