#!/usr/bin/perl
use strict;
use warnings;

=head1 NAME

bp_translate_feature - translates features

=head1 SYNOPSIS

bp_translate_feature E<lt> cdna_cds.fa E<gt> protein.fa

=head1 DESCRIPTION 

The script will translate one fasta file (on stdin) to protein on stdout

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
email or the web:

  https://redmine.open-bio.org/projects/bioperl/

=head1 AUTHOR

  Ewan Birney E<lt>birney@ebi.ac.ukE<gt>

=cut

use Bio::FeatureIO;
use Getopt::Long;

my ($format) = 'gff';

GetOptions(
	   'format:s'  => \$format,
	   );

my $oformat = 'faldo';

# this implicity uses the <> file stream
my $featurein = Bio::FeatureIO->new( -format => $format, -file => shift); 
my $featureout = Bio::FeatureIO->new( -format => $oformat, -file => ">-" );


while ( my @f = $featurein->next_feature_group ) {
    $featureout->write_feature($_) for @f;
}

__END__
