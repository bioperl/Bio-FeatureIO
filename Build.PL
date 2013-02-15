#!/usr/bin/perl -w

use strict;
use warnings;
use Module::Build;

my $build = Module::Build->new(
    dist_name            => 'Bio::FeatureIO',
    dist_version         => '0.01',
    module_name          => 'Bio::FeatureIO',
    dist_author          => 'BioPerl Team <bioperl-l@bioperl.org>',
    dist_abstract        => 'Bio::FeatureIO - modules for reading and writing sequence features',
    license              => 'perl',
    recursive_test_files => 1,
    create_license       => 1,
    meta_merge           => {
	                     'resources' => {
				         'repository' => 'http://github.com/bioperl/Bio-FeatureIO',
                             },
			    },
    auto_features       => {
                            Network => {
                                        description => "Enable tests that need an internet connection",
                                        requires    => { 'LWP::UserAgent' => 0 }
                                       }
                           },
    get_options         => {
                            accept  => { },
                            network => { }
                           },
    requires             => {
                             'perl'                => '5.10.0',
                             'Bio::Root::Version'  => '1.006900',
                             'Bio::GFF3::LowLevel' => 0,
                             'IO::String'          => 0,
                             'URI::Escape'         => '0',    # dealing with web resources/Bio::FeatureIO::gff,Bio::FeatureIO::interpro,Bio::DB::Biblio::eutils,Bio::DB::EUtilParameters,Bio::DB::Query::GenBank,Bio::DB::NCBIHelper,Bio::SeqFeature::Annotated',
                             'XML::DOM::XPath'     => '0.13', # parsing interpro features/Bio::FeatureIO::interpro',
                            },
    build_requires       => {'Test::Warn' => 0 },
    recursive_test_files => 1,
    dynamic_config       => 1,
);

my $accept = $build->args->{accept};

# Do network tests?
my $do_network_tests = 0;
if ($build->args('network')) {
    $do_network_tests = $build->feature('Network');
}
elsif ($build->feature('Network')) {
    $do_network_tests = $accept ? 0 :
    $build->y_n("Do you want to run tests that require connection to ".
        "servers across the internet\n(likely to cause some failures)? y/n",
        'n');
}

if ($do_network_tests) {
    $build->notes(network => 1);
    $build->log_info("  - will run internet-requiring tests\n");
} else {
    $build->notes(network => 0);
    $build->log_info("  - will not run internet-requiring tests\n");
}

# Create the build script and exit
$build->create_build_script;

exit;
