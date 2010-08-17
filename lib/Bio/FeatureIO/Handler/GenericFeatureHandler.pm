package Bio::FeatureIO::Handler::GenericFeatureHandler;

use base qw(Bio::Root::Root Bio::HandlerBaseI);

use strict;
use warnings;
use Data::Dumper;
use Bio::SeqFeature::Generic;
use Bio::SeqFeature::Tools::Unflattener;
use Bio::SeqIO;

my $ct = 0;
my %GFF3_RESERVED_TAGS = map {$_ => $ct++ }
    qw(ID Name Alias Parent Target Gap
    Derives_from Note Dbxref Ontology_term Index);
    
my %HANDLERS = (
    'directive'             => \&directives,
    'comment'               => \&comment,
    'feature'               => \&seqfeature,
    'sequence'              => \&sequence,
);

our $ONTOLOGY_STORE;
our $UNFLATTENER;
our $ID_HANDLER;

sub new {
    my ($class, @args) = @_;
    my $self = $class->SUPER::new(@args);
    $self = {@args};
    bless $self,$class;
    $self->handler_methods();
    return $self;
}

sub data_handler {
    my ($self, $data) = @_;
    my $nm = $data->{MODE} || $self->throw("No type tag defined!\n".Dumper($data));
    
    $self->set_parameters('mode', $nm eq 'directive' ? $data->{DATA}->{type} : $nm);
    
    my $method = (exists $self->{'handlers'}->{$nm}) ? ($self->{'handlers'}->{$nm}) :
                (exists $self->{'handlers'}->{'_DEFAULT_'}) ? ($self->{'handlers'}->{'_DEFAULT_'}) :
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

sub set_handler_helper {
    my ($self, $name, $sub) = @_;
    return if !($name && $sub);
    $self->throw("Passed callback must be a code ref") if $sub && ref $sub eq 'CODE';
    $self->{'handlers'}->{$name} = $sub;
}

sub format {
    my $self = shift;
    return $self->{format} = shift if @_;
    return $self->{format};
}

sub fast {
    my $self = shift;
    $self->{fast} || 0;
}

sub reset_parameters {
    my ($self) = @_;
    $self->{parameters} = {};
}

sub get_parameters {
    my ($self, $param) = @_;
    return if !($param);
    $self->{parameters}->{$param};
}

sub set_parameters {
    my ($self, $param, $value) = @_;
    return if !($param && defined $value);
    $self->{parameters}->{$param} = $value;
}

# this needs to be a Bio::SeqFeature::CollectionI that can distinguish
# between sequence regions; the simplest versions don't

sub feature_collection {
    my $self = shift;
    return $self->{feature_collection} = shift if @_;
    return $self->{feature_collection};
}

sub file_handle {
    return shift->{-fh};
}

sub fasta_mode {
    my $self = shift;
    my $mode = $self->get_parameters('mode');
    return unless $mode;
    $mode eq 'sequence' || $mode eq 'sequence-region' ? 1 : 0;
}

sub resolve_references {
    my $self = shift;
    my $mode = $self->get_parameters('mode');
    return unless $mode;
    $mode eq 'resolve-references' || $mode eq 'sequence' || $mode eq 'sequence-region' ? 1 : 0;
}

################ HANDLERS ################

# Handler methods are designed so they are called as sub references, not as
# class or instance based calls. This decouples them from any handler class and
# allow some customization (for instance, if we need the parent parser to
# override them). The parser in this case is ultimately responsible for
# determining what happens to the data.

# Note this just passes in the data w/o munging it beyond recognition
sub seqfeature {
    my ($handler, $data) = @_;

    my %sf_data = map {'-'.$_ => $data->{DATA}->{$_}}
        grep { $data->{DATA}->{$_} ne '.' }
        sort keys %{$data->{DATA}};
    
    if ($data->{DATA}->{attributes}) {
        delete $sf_data{-attributes};
        my %tags;
        
        # TODO: GFF3-specific split; need to make more general
        for my $kv (split(/\s*;\s*/, $data->{DATA}->{attributes})) {
            my ($key, $rest) = split(/[=\s]/, $kv, 2);
            # add optional/required URI unescaping here
            my @vals = split(',',$rest);
            $tags{$key} = \@vals;
        }
        $sf_data{-tag} = \%tags;
    }
    
    return Bio::SeqFeature::Generic->new(%sf_data);
}

sub directives {
    my ($handler, $data) = @_;
    my $directive = $data->{DATA}->{type};
    if ($directive eq 'sequence') {
        my $fh = $handler->file_handle;
        $handler->throw("Handler doesn't have a set file handle") if !$fh;
        return Bio::SeqIO->new(-format => 'fasta',
                       -fh     => $fh);
    } elsif ($directive eq 'sequence-region') {
        # we can make returning a feature optional here, but we should do
        # something with the data in all cases
        
        my $sf_data = $data->{DATA};
        return Bio::SeqFeature::Generic->new(-start     => $sf_data->{start},
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
    # if we reach this point, the sequence stream has already been read, so
    # we need to seek back to the start point.  Note if the stream isn't seekable
    # this will fail spectacularly at this point!
    my ($start, $len) = @{$data}{qw(START LENGTH)};
    my $fh = $handler->file_handle;
    $handler->throw("Handler doesn't have a set file handle") if !$fh;
    seek($fh, $start, 0);
    return Bio::SeqIO->new(-format => 'fasta',
                           -fh     => $fh);
}

# no op, we just skip these
sub comment {}

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
  
L<bioperl-l@bioperl.org>
  
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
