# -*-Perl-*- Test Harness script for Bioperl
# $Id: FeatureIO.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;

BEGIN {
    use lib './inc';
    use Bio::Root::Test;

    test_begin(-tests => 4);

    use_ok($_) for qw(
        Bio::FeatureIO
        Bio::FeatureIO::gff
        Bio::FeatureIO::ptt
        Bio::FeatureIO::vecscreen_simple
    );
}

done_testing();

exit;
