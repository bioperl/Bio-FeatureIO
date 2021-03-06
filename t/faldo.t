# -*-Perl-*- Test Harness script for Bioperl
# $Id: FeatureIO.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib './inc';
    use Bio::Root::Test;
    use Bio::FeatureIO;
    use Bio::GFF3::LowLevel 'gff3_parse_feature';
    use IO::String;

    test_begin(-tests => 4);
}

my ($io, $f, $s, $fcount, $scount);

################################################################################
#
# translate from GFF3 to FALDO
#

for my $modifier ( sub {}, sub { $_[0]->remove_tag('Parent') if $_[0]->has_tag('Parent') } ) {
    my $kg_file = test_input_file('tomato_test.gff3');
    ok( $io = Bio::FeatureIO->new( -file => $kg_file ) );
    my $out_faldo = '';
    my $out = Bio::FeatureIO->new(
        -fh      => IO::String->new( \$out_faldo ),
        -format  => 'faldo',
        -version => 3,
      );

    # we just 'trust' the output for now - we have no faldo reader
    ok( $out );
    while ( my @f = $io->next_feature_group ) {
        #$modifier->($_) for @f, map $_->get_SeqFeatures, @f;
        $out->write_feature($_) for @f;
    }
    print $out_faldo;


}
done_testing();
exit;
