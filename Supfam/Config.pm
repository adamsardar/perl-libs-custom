package Supfam::Config;
require Exporter;

our @ISA       = qw(Exporter);
our @EXPORT    = qw(%SUPFAM %CONFIG);
our @EXPORT_OK = qw();
our $VERSION   = 1.00;

use strict;
use warnings;
#use diagnostics;
use XML::Simple qw(:strict); #Used to load in the config file and make a hash.

eval {
   our %SUPFAM = %{XMLin($ENV{'HOME'}."/supfam_config.xml", ForceArray => 0, KeyAttr => [ ])};
   our %CONFIG = %{XMLin("config.xml", ForceArray => 0, KeyAttr => [ ])};
};

#Exception handling stub for when we dont have a config file
if ($@) {
   die "Fatal Error: No supfam_config.xml found. $@\n" if $@ =~ /supfam_config.xml/;
   #Leave it to plugin modules to do: if (not defined $CONFIG) but still warn of possible problems
  # warn "Not using a local config.xml for this script.\n" if ($@ =~ /config.xml/);
}

1;
__END__

=head1 NAME

Supfam::Config.pm

=head1 DESCRIPTION

Provides configuration information for the SUPERFAMILY database and related databases. Loads in the data from lib/Supfam/supfam_config.xml as well as any local conf.xml
in the working directory of the currently executing script.

=head2 $SUPFAM

Pointer to a configuration hash with similar structure to the XML definition found in the supfam_config.xml file.
Example use: `$SUPFAM->{'database'}'name'}`

=head2 $CONFIG

Pointer to a configuration hash with similar structure to the XML definition found in the local config.xml file.

=cut
