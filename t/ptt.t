# -*-Perl-*- Test Harness script for Bioperl
# $Id: FeatureIO.t 15112 2008-12-08 18:12:38Z sendu $

use strict;
use warnings;
use lib './inc';
use Bio::Root::Test;
use Bio::FeatureIO;
use Data::Dumper;

my ($io, $f, $s, $fcount, $scount);

$fcount = 0;

my $ptt_in = Bio::FeatureIO->new(
  -file => test_input_file('test.ptt'), 
  -format => 'ptt',
);
ok($ptt_in);

while (my $f = $ptt_in->next_feature) {
    $fcount++;
    if ($fcount==2) {
        # 2491..3423  + 310 24217063  metF  LB002 - COG0685E  5,10-methylenetetrahydrofolate reductase
        is( $f->start , 2491 );
        is( $f->end , 3423 );
        is( $f->strand, 1);
        is( $f->display_name, 'metF');
        is( $f->primary_id, '24217063');
        is( ($f->get_tag_values('ID'))[0],'24217063' );
        is( ($f->get_tag_values('Alias'))[0], 'LB002' );
        ok( ! $f->has_tag('code') );
        is( ($f->get_tag_values('cog'))[0],'COG0685E' );
        is( ($f->get_tag_values('product'))[0], '5,10-methylenetetrahydrofolate reductase' );   
    }
}
is($fcount, 367, 'ptt file count');

done_testing();

exit;
