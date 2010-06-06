# $Id: SimpleEventGeneratorI.pm 16108 2009-09-16 17:07:49Z cjfields $
#
# BioPerl module for Bio::Event::SimpleEventGeneratorI
#
# Please direct questions and support issues to <bioperl-l@bioperl.org> 
#
# Cared for by Chris Fields <cjfields at bioperl dot org>
#
# Copyright Chris Fields
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::Event::SimpleEventParserI- This interface describes a parser for generating
a simple stream of events.  

=head1 SYNOPSIS

    # Do not use this object directly
    # This object has the basic methods for describing an event generator

=head1 DESCRIPTION

This object augments the basic EventGeneratorI interface with additional methods
that indicate a stream of simple, non-object-based events is generated from a
raw data stream. These are generally in the form of hash references containing
one or more defined types of data and can be as simple as single pieces of data
(analagous to a SAX-based parser) or contain pre-processed chunks of data that
can either be called on directly by the user for faster processing or be passed
into a common Bio::Event::EventHandlerI for further processing into various
BioPerl objects.

The main stream interface (for example, Bio::FeatureIO) or plugins of a
particular stream interface that use this interface (Bio::SeqIO::foo classes)
should define the data types using a common set of terms, possibly defined in an
ontology; similarly, the event handler should defined what data is acceptable.

In particular, the following formats are suggested:
    
   # SAX-like, no additional kv pairs, can be empty DATA hash
   $VAR = {
       'MODE'      => 'FOO',
       'DATA'      => {
           'DATA  => 'FOO', 
       }
   ;
    
or
    
   # processed, kv pairs
   $VAR = {
       'MODE'      => 'FOO',
       'DATA'      => {
           # simple 
           'KEY1' => 'VAL',
           
           # multiple values per key
           'KEY2' => ['VAL1', 'VAL2'],  # multiples
           
           # hierarchy
           'KEY3' => {'TAG1' => 'VAL'},
           'KEY4' => {'TAG2' => ['VAL1','VAL2'] },
       }
   };

For convenience, in the second version above key names can be parameters
passed directly into a new instance; this will be largely dictated by the
requirements of the event handler, but allows for some flexibility when
customizing data:

   while (my $data = $iterator->next_dataset) {
      # do some damage to $data
    
      my $sf = Bio::SeqFeature::Foo->new(%$data);
   }

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Chris Fields

Email cjfields at bioperl dot org

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...

package Bio::Event::SimpleEventGeneratorI;
use strict;

use base qw(Bio::Event::EventGeneratorI);

sub next_dataset {
    shift->throw_not_implemented;
}

sub handler {
    shift->throw_not_implemented;
}

1;
