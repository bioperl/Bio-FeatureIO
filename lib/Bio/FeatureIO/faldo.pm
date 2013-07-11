package Bio::FeatureIO::faldo;

use strict;
use warnings;
use 5.010;
use base qw(Bio::FeatureIO);
use Bio::FeatureIO::SequenceOntologyTypeMap;

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    $self->{_typemap} = Bio::FeatureIO::SequenceOntologyTypeMap->get_typemap();
    $self->{_prefixmap_written} = {};
    

}


=head2 write_feature()

 Usage   : $featureio->write_feature( Bio::SeqFeature::Annotated->new(...) );
 Function: writes a feature in FALDO format.  the FALDO version used is
           governed by the '-version' argument passed to Bio::FeatureIO->new(),
           and defaults to FALDO version 3.
 Returns : nothing meaningful
 Args    : a Bio::SeqFeature::Annotated object.

=cut

sub write_feature {
    my($self,$feature) = @_;

    unless( $self->{wrote_ttl_header} ) {
        $self->write_ttl_header();
        $self->{wrote_ttl_header} = 1;
    }

    unless( $feature ) {
        $self->throw("faldo.pm cannot write_feature unless you give a feature to write.\n");
    }


    $self->_write_feature_with_parent( $feature );
    $self->_write_feature_with_parent( $_, $feature ) for $feature->get_SeqFeatures;


}

sub _write_feature_with_parent {
    my($self,$f, $parent) = @_;
    
    # type
    my ($id) = $f->get_tag_values('ID');
    my $uri = $self->_id_to_uri($id);
    $self->_write_triple($uri, 'rdf:type', $self->_so($f->primary_tag));

    # begin and end
    my $begin_obj = $self->_skolem($uri, 'b');
    my $end_obj = $self->_skolem($uri, 'e');
    $self->_write_triple($begin_obj, 'rdf:type', 'faldo:ExactPosition');
    $self->_write_triple($begin_obj, 'faldo:position', _xsd( int=> $f->start ));
    $self->_write_triple($end_obj, 'rdf:type', 'faldo:ExactPosition');
    $self->_write_triple($end_obj, 'faldo:position', _xsd( int=> $f->end ));
    $self->_write_triple($uri, 'faldo:begin', $begin_obj);
    $self->_write_triple($uri, 'faldo:end', $end_obj);

    if ($parent) {
        # TODO
    }

    $self->_print("\n");
    
    
}

sub _so {
    my ($self, $t) = @_;

    if ($self->{_prefixmap_written}->{$t}) {
        return "$t:";
    }
    my $so_id = $self->{_typemap}->{$t};
    if ($so_id) {
        $so_id =~ s/^SO:/obo:SO_/;
        return $so_id
    }
    else {
        $self->warn("No typemap for "+$t);
        return "obo:SO_".$t;
    }
}

sub uri_prefix_for_features {
   my ($self, $id) = @_;
   # TODO: configurability
   return "http://bioperl.org/features/";
}

sub _id_to_uri {
    my ($self, $id) = @_;
    my $pfx = $self->uri_prefix_for_features($id);
    # TODO: configurability
    return "<$pfx$id>";
}

sub _skolem {
    my ($self, $base, $suffix) = @_;
    if ($base =~ m@^\<(.*)\>@) {
        return "<$1"."__$suffix>";
    }
    else {
        return "$base"."__$suffix";
    }
}

sub _xsd {
    my ($t,$v) = @_;
    $v =~ s/\"/\\\"/g; # TODO - check this against TTL grammar
    return '"'.$v.'"^^xsd:'.$t;
}

# we assume the arguments have already been translated to ttl
# syntax for URIs or literals
sub _write_triple {
    my ($self,$s,$p,$o) = @_;
    $self->_print( "$s $p $o .\n");
    return;
}

sub _uri {
    my ($self,$obj) = @_;
    return "$obj";
}

sub _prefixmap {
    return 
      (rdf => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
       rdfs => 'http://www.w3.org/2000/01/rdf-schema#',
       xsd => 'http://www.w3.org/2001/XMLSchema#',
       obo => 'http://purl.obolibrary.org/obo/',
       faldo => 'http://biohackathon.org/resource/faldo#',
      );

}

sub write_ttl_header {
    my ($self) = @_;
    my %pm = $self->_prefixmap;
    
    foreach my $k (keys %pm) {
        $self->_write_prefix_declaration($k, $pm{$k});
    }

    # use shorthand for features
    foreach my $t ($self->get_feature_types_used) {
        $self->_write_prefix_declaration($t, $self->_so($t));
    }

    $self->_print("\n");
}

sub get_feature_types_used {
    my ($self) = @_;
    return qw(gene transcript exon CDS intron mRNA five_prime_UTR three_prime_UTR);
}

sub _write_prefix_declaration {
    my ($self, $k, $v) = @_;
    $self->_print('@prefix '.$k.': <'.$v.'> .'."\n");
    $self->{_prefixmap_written}->{$k} = $v;
}

sub _write_feature_3 {
    my ( $self, $feature, $parent_feature ) = @_;


}

sub _faldo3_lowlevel_hashref {
    my ( $self, $f, $parent ) = @_;

    my @tags = $f->get_all_tags;
    if( $f->can('phase') ) {
        @tags = grep $_ ne 'phase', @tags;
    }

    if( $parent && ! $parent->has_tag('ID') ) {
        $self->throw("feature ".($parent->display_name||'(no display name)')." has subfeatures but no 'ID' tag, cannot write as FALDO3");
    }
    my $parent_id = $parent ? ($parent->get_tag_values('ID'))[0] : ();

    return {
        seq_id => $f->seq_id,
        source => $f->source_tag,
        type   => $f->primary_tag,
        start  => $f->start,
        end    => $f->end,
        score  => $f->score,
        strand => ( $f->strand ? $f->strand == 1 ? '+' : '-'
                                : undef
                  ),
        phase  => ($f->can('phase') ? $f->phase : undef ),
        attributes => {
            ( map {
                my $tag = $_;
                $tag => [ $f->get_tag_values( $tag ) ]
              } @tags,
            ),
            ( $parent_id ? ( Parent => [ $parent_id ] ) : () ),
          },
    };
}

1;

__END__

=head1 NAME

Bio::FeatureIO::faldo - read/write FALDO feature files

=head1 SYNOPSIS

  my $feature; #get a Bio::SeqFeatureI somehow
  my $featureOut = Bio::FeatureIO->new(
    -format => 'faldo',
    -fh => \*STDOUT,
  );
  $featureOut->write_feature($feature);

=head1 DESCRIPTION

This is a write for FALDO. Parsing is not yet implemented.


=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org                 - General discussion
  http://bioperl.org/wiki/Mailing_list  - About the mailing lists

=head2 Support 

Please direct usage questions or support issues to the mailing list:

I<bioperl-l@bioperl.org>

rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.open-bio.org/

=head1 AUTHORS

 Chris Fields, <cjfields at bioperl dot org>
 Robert Buels, <rbuels at cpan dot org>

Refactored from the original work by:

 Allen Day, <allenday@ucla.edu>

=head1 CONTRIBUTORS

 Steffen Grossmann, <grossman@molgen.mpg.de>
 Scott Cain, <scain@cpan.org>
 Rob Edwards <rob@salmonella.org>

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=head2 next_feature()

 Usage   : my $feature = $featureio->next_feature();
 Function: reads a feature record from a FALDO stream and returns it as an object.
 Returns : a Bio::SeqFeature::Annotated object
 Args    : N/A

=head2 next_feature_group

 Usage   : @feature_group = $stream->next_feature_group
 Function: Reads the next feature_group from $stream and returns it.

           Feature groups in FALDO3 files are separated by '###' directives. The
           features in a group might form a hierarchical structure. The
           complete hierarchy of features is returned, i.e. the returned array
           represents only the top-level features.  Lower-level features can
           be accessed using the 'get_SeqFeatures' method recursively.

 Example : # getting the complete hierarchy of features in a FALDO3 file
           my @toplevel_features;
           while (my @fg = $stream->next_feature_group) {
               push(@toplevel_features, @fg);
           }
 Returns : an array of Bio::SeqFeature::Annotated objects
 Args    : none

=head2 next_seq()

  Usage   : $featureio->next_seq( );
  Function: access the FASTA section (if any) at the end of the FALDO stream. Note
            that this method will return undef before all the features in the FALDO
            stream have been handled.
  Returns : a Bio::SeqI object or undef
  Args    : none

=cut
    
