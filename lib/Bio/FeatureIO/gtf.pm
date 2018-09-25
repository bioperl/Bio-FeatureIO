package Bio::FeatureIO::gtf;

use utf8;
use strict;
use warnings;
use base qw(Bio::FeatureIO::gff);

# ABSTRACT: read write features in GTF format
# AUTHOR:   Allen Day <allenday@ucla.edu>
# OWNER:    Allen Day
# LICENSE:  Perl_5

# Object preamble - inherits from Bio::Root::Root

sub _initialize {
    my($self,%arg) = @_;
    $arg{-version} = 2.5;
    $self->SUPER::_initialize(%arg);
    return 1;
}



1;

__END__

=head1 SYNOPSIS

Please see L<Bio::FeatureIO::gff>

=head1 DESCRIPTION

GTF, is also known as GFF v2.5.  This class is simply a subclass
of Bio::FeatureIO::gff that initializes with -version =E<gt> 2.5.

=cut
