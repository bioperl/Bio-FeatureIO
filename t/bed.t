# -*-Perl-*- Test Harness script for Bioperl

use strict;
use warnings;
use lib './inc';
use Bio::Root::Test;
use Bio::FeatureIO;
use Data::Dumper;

my ($io, $f, $s, $fcount, $scount);

################################################################################
#
# use FeatureIO::bed to read a bed file
#

ok($io = Bio::FeatureIO->new(-file => test_input_file('hg19_sample.bed')));

ok($f = $io->next_feature);

# Check correct conversion of [0, feature-end+1) bed-coordinates into [1, feature-end]
# bioperl coordinates.  (here: bed [0, 10))
is($f->start, 23281095);
is($f->end, 23283550);
is($f->strand, 1);

# Check field values.
is($f->display_name, "BC043238");
is($f->seq_id, "chr2");

#TODO: blocks into subfeatures

ok($f = $io->next_feature);
is($f->start, 23348659);
is($f->end, 23349587);
is($f->strand, -1);

# Check field values.
is($f->display_name, "AK024624");
is($f->seq_id, "chr2");

my $ct;

# get rest of features
while ($io->next_feature) {
    $ct++;
}

is($ct, 6);

ok($io = Bio::FeatureIO->new(-file => test_input_file('1.bed')));
ok($f = $io->next_feature);

is($f->start, 1);
is($f->end, 10);

# Check field values.
is($f->display_name, "test-coordinates-1");
is($f->seq_id, "chr1");

done_testing();

exit;
