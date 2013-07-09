package Bio::FeatureIO::Handler::GenericFeatureHandler;

use base qw(Bio::Root::Root Bio::HandlerBaseI);

use 5.010;
use strict;
use warnings;
use Bio::Factory::ObjectFactory;
use Bio::SeqIO;
use Data::Dumper;

my $ct = 0;
my %GFF3_RESERVED_TAGS = map {$_ => $ct++ }
    qw(ID Name Alias Parent Target Gap
    Derives_from Note Dbxref Ontology_term Index);
    
my %HANDLERS = (
    # GFF3-specific
    'directive'             => \&directives,
    
    # BED-specific
    'track_definition'      => \&track_definition,
    
    # generic
    'comment'               => \&comment,
    'feature'               => \&seqfeature,
    'sequence'              => \&sequence,
    
    # defaults (none for now)
    #'default'               => \&default
);

our $ONTOLOGY_STORE;
our $UNFLATTENER;
our $ID_HANDLER;

sub new {
    my ($class, %args) = @_;
    my %atts = map {(my $newk = $_ ) =~ s/^-//; $newk => $args{$_} } keys %args;
    my $self = bless \%atts, $class;
    return $self;
}

sub data_handler {
    my ($self, $data) = @_;
    my $nm = $data->{MODE} || $self->throw("No type tag defined!");

    $self->{current_state} = $data;

    my $method = (exists $HANDLERS{$nm}) ? ($HANDLERS{$nm}) :
                (exists $HANDLERS{'_DEFAULT_'}) ? ($HANDLERS{'_DEFAULT_'}) :
                undef;
    
    if ($method && ref $method eq 'CODE') {
        return $self->$method($data);
    } else {
        $self->debug("No handler defined for $nm\n");
        return;
    }
}

sub handler_methods {
    my $self = shift;
    my $handlers = shift;
    if (!($self->{'handlers'}) || defined $handlers && ref $handlers eq 'HASH') {
        $self->{'handlers'} = $handlers || \%HANDLERS;
    }
    return ($self->{'handlers'});
}

sub format {
    my $self = shift;
    return $self->{-format} = shift if @_;
    return $self->{-format};
}

sub fast {
    my $self = shift;
    $self->{-fast} || 0;
}

sub get_parameter {
    my ($self, $param) = @_;
    return if !($param);
    $self->{parameters}{$param};
}

*get_parameters = \&get_parameter;

sub set_parameter {
    my ($self, $param, $value) = @_;
    return if !($param && defined $value);
    $self->{parameters}{$param} = $value;
}

*set_parameters = \&set_parameter;

sub get_parameter_names {
    my $self = shift;
    return sort keys %{$self->{parameters}} if $self->{parameters};
    return ();
}

sub reset_parameters {
    my ($self, $param) = @_;
    $self->{parameters} = {};
}

sub file_handle {
    my ($self, $fh) = @_;
    return $self->{-fh} = $fh if $fh;
    return $self->{-fh};
}

# GFF3-specific convenience methods
sub fasta_mode {
    my $self = shift;
    return unless $self->{current_state}; # init
    my ($mode, $dir_type) = ($self->{current_state}{MODE} || '',
                             $self->{current_state}{DATA}{type} || '');
    $mode eq 'sequence' || $dir_type eq 'sequence';
}

sub resolve_references {
    my $self = shift;
    my $dir_type = $self->{current_state}{DATA}{type} || '';
    $dir_type eq 'resolve-references' || $self->fasta_mode;
}

sub feature_factory {
    my ($self, $factory) = @_;
    if ($factory && $factory->isa('Bio::Factory::ObjectFactoryI')) {
        $self->{feature_factory} = $factory;
    }
    if (! defined $self->{feature_factory}) {
        $self->{feature_factory} = Bio::Factory::ObjectFactory->new(
                                -type       => 'Bio::SeqFeature::Generic',
                                -interface  => 'Bio::SeqFeatureI'    );
    }
    $self->{feature_factory}
}

sub feature_class {
    my $self = shift;
    return $self->{-feature_class} = shift if @_;
    $self->{-feature_class} || 'Bio::SeqFeature::Generic';
}

################ HANDLERS ################

# Handler methods are designed so they are called as sub references, not as
# class or instance based calls. This decouples them from any handler class and
# allow some customization (for instance, if we need the parent parser to
# override them). The parser in this case is ultimately responsible for
# determining what happens to the data.

sub seqfeature {
    my ($handler, $data) = @_;
    # TODO: Optionally type check ontology here, preferably prior to expending
    #       unnecessary energy/time on creating objects (Robert Bradbury rules)
    return $handler->feature_factory->create_object(%{$data->{DATA}});
}

sub directives {
    my ($handler, $data) = @_;
    my $directive = $data->{DATA}->{type};
    return unless $directive;
    if ($directive eq 'sequence') {
        my $fh = $handler->file_handle;
        $handler->throw("Handler doesn't have a set file handle") if !$fh;
        return Bio::SeqIO->new(-format => 'fasta',
                       -fh     => $fh);
    } elsif ($directive eq 'sequence-region') {
        # we can make returning a feature optional here, but we should do
        # something with the data in all cases
        my $sf_data = $data->{DATA};
        return $handler->feature_factory->create_object(
            -start     => $sf_data->{start},
            -end       => $sf_data->{end},
            -strand    => 1,
            -seq_id    => $sf_data->{id},
            -primary_tag  => 'region');
    } else {
        # defaults for other directives?
    }
    return;
}

sub sequence {
    my ($handler, $data) = @_;
    
    # So, at this point we have to decide whether to indicate this is the start
    # of a sequence stream (i.e. grab the file handle and create a SeqIO), or
    # allow further sequence data to be passed in. We punt and do the easy thing
    # for now, but we should allow flexibility here, so maybe something
    # configurable?
    
    # Note this relies on having knowledge of the specific place in the stream
    
    my ($start, $len) = @{$data}{qw(START LENGTH)};
    my $fh = $handler->file_handle;
    $handler->throw("Handler doesn't have a set file handle") if !$fh;
    seek($fh, $start, 0);
    return Bio::SeqIO->new(-format => 'fasta',
                           -fh     => $fh);
}

# no ops, we just skip these
sub comment {}

sub track_definition {}

1;

__END__

=head1 NAME

Bio::FeatureIO::Handler::GenericFeatureHandler.pm - <One-line description of module's
purpose>

=head1 VERSION

This documentation refers to Bio::FeatureIO::Handler::GenericFeatureHandler.pm version
Bio::Root::Root.

=head1 SYNOPSIS

   use Bio::FeatureIO::Handler::GenericFeatureHandler.pm;
   # Brief but working code example(s) here showing the most common usage(s)

   # This section will be as far as many users bother reading,

   # so make it as educational and exemplary as possible.

=head1 DESCRIPTION

<TODO>
A full description of the module and its features.
May include numerous subsections (i.e., =head2, =head3, etc.).

=head1 SUBROUTINES/METHODS

<TODO>
A separate section listing the public components of the module's interface.
These normally consist of either subroutines that may be exported, or methods
that may be called on objects belonging to the classes that the module provides.
Name the section accordingly.

In an object-oriented module, this section should begin with a sentence of the
form "An object of this class represents...", to give the reader a high-level
context to help them understand the methods that are subsequently described.

=head1 DIAGNOSTICS

<TODO>
A list of every error and warning message that the module can generate
(even the ones that will "never happen"), with a full explanation of each
problem, one or more likely causes, and any suggested remedies.

=head1 CONFIGURATION AND ENVIRONMENT

<TODO>
A full explanation of any configuration system(s) used by the module,
including the names and locations of any configuration files, and the
meaning of any environment variables or properties that can be set. These
descriptions must also include details of any configuration language used.

=head1 DEPENDENCIES

<TODO>
A list of all the other modules that this module relies upon, including any
restrictions on versions, and an indication of whether these required modules are
part of the standard Perl distribution, part of the module's distribution,
or must be installed separately.

=head1 INCOMPATIBILITIES

<TODO>
A list of any modules that this module cannot be used in conjunction with.
This may be due to name conflicts in the interface, or competition for
system or program resources, or due to internal limitations of Perl
(for example, many modules that use source code filters are mutually
incompatible).

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

User feedback is an integral part of the evolution of this and other Biome and
BioPerl modules. Send your comments and suggestions preferably to one of the
BioPerl mailing lists. Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

Patches are always welcome.

=head2 Support 

Please direct usage questions or support issues to the mailing list:

 bioperl-l@bioperl.org

rather than to the module maintainer directly. Many experienced and reponsive
experts will be able look at the problem and quickly address it. Please include
a thorough description of the problem with code and data examples if at all
possible.

=head2 Reporting Bugs

Preferrably, Biome bug reports should be reported to the GitHub Issues bug
tracking system:

  http://github.com/cjfields/biome/issues

Bugs can also be reported using the BioPerl bug tracking system, submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 EXAMPLES

<TODO>
Many people learn better by example than by explanation, and most learn better
by a combination of the two. Providing a /demo directory stocked with
well-commented examples is an excellent idea, but your users might not have
access to the original distribution, and the demos are unlikely to have been
installed for them. Adding a few illustrative examples in the documentation
itself can greatly increase the "learnability" of your code.

=head1 FREQUENTLY ASKED QUESTIONS

<TODO>
Incorporating a list of correct answers to common questions may seem like extra
work (especially when it comes to maintaining that list), but in many cases it
actually saves time. Frequently asked questions are frequently emailed
questions, and you already have too much email to deal with. If you find
yourself repeatedly answering the same question by email, in a newsgroup, on a
web site, or in person, answer that question in your documentation as well. Not
only is this likely to reduce the number of queries on that topic you
subsequently receive, it also means that anyone who does ask you directly can
simply be directed to read the fine manual.

=head1 COMMON USAGE MISTAKES

<TODO>
This section is really "Frequently Unasked Questions". With just about any kind
of software, people inevitably misunderstand the same concepts and misuse the
same components. By drawing attention to these common errors, explaining the
misconceptions involved, and pointing out the correct alternatives, you can once
again pre-empt a large amount of unproductive correspondence. Perl itself
provides documentation of this kind, in the form of the perltrap manpage.

=head1 SEE ALSO

<TODO>
Often there will be other modules and applications that are possible
alternatives to using your software. Or other documentation that would be of use
to the users of your software. Or a journal article or book that explains the
ideas on which the software is based. Listing those in a "See Also" section
allows people to understand your software better and to find the best solution
for their problem themselves, without asking you directly.

By now you have no doubt detected the ulterior motive for providing more
extensive user manuals and written advice. User documentation is all about not
having to actually talk to users.

=head1 (DISCLAIMER OF) WARRANTY

<TODO>
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=head1 ACKNOWLEDGEMENTS

<TODO>
Acknowledging any help you received in developing and improving your software is
plain good manners. But expressing your appreciation isn't only courteous; it's
also enlightened self-interest. Inevitably people will send you bug reports for
your software. But what you'd much prefer them to send you are bug reports
accompanied by working bug fixes. Publicly thanking those who have already done
that in the past is a great way to remind people that patches are always
welcome.

=head1 AUTHOR

Chris Fields  (cjfields at bioperl dot org)

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2010 Chris Fields (cjfields at bioperl dot org). All rights reserved.

followed by whatever licence you wish to release it under.
For Perl code that is often just:

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.
