use strict;

use Test::More tests => 5;

use_ok('Bio::SeqFeature::Annotated');

my $fea = Bio::SeqFeature::Annotated->new();
isa_ok($fea, "Bio::SeqFeatureI",'isa SeqFeatureI');
isa_ok($fea, "Bio::AnnotatableI",'isa AnnotatableI');
$fea = Bio::SeqFeature::Generic->new();
isa_ok($fea, "Bio::SeqFeatureI",'isa SeqFeatureI');
isa_ok($fea, "Bio::AnnotatableI",'isa AnnotatableI');
