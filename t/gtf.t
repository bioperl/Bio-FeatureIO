# -*-Perl-*- Test Harness script for Bioperl
# $Id: FeatureIO.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;
use lib './inc';
use Bio::Root::Test;
use Bio::FeatureIO;

# this is mainly GFF3-specific, GFF2/GTF to be added

my ($io, $f, $s, $fcount, $scount);

################################################################################
#
# use FeatureIO::gff to read a FASTA file.
#

$fcount = 0;
$scount = 0;

ok( $io = Bio::FeatureIO->new(-format => 'gtf',
                              -file => test_input_file('dna1.fa') ) );

is ($io->version, '2.5');

#read features
while($f = $io->next_feature()){
    $fcount++;
}
is($fcount, 0);

# GTF do not have sequences associated with them, but if a sequence is provided
# it will read it
while($s = $io->next_seq()){
    $scount++;
    if ($scount == 1) {
        is($s->id, 'Test1');
    }
}
is($scount,  1);

ok( $io = Bio::FeatureIO->new(-format => 'gtf',
                              -file => test_input_file('mm9_sample_ucsc.gtf') ) );

while($f = $io->next_feature()){
    $fcount++;
}
is($fcount, 55);

done_testing();

