package Bio::FeatureIO::faldo;

use strict;
use warnings;
use 5.010;
use base qw(Bio::FeatureIO);
use Bio::FeatureIO::OntologyMapper;

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    $self->{_typemap} = Bio::FeatureIO::OntologyMapper->get_typemap();
    $self->{_prefixmap_written} = {};
    

}

sub _prefixmap {
    return 
      (rdf => 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
       rdfs => 'http://www.w3.org/2000/01/rdf-schema#',
       xsd => 'http://www.w3.org/2001/XMLSchema#',
       obo => 'http://purl.obolibrary.org/obo/',
       faldo => 'http://biohackathon.org/resource/faldo#',
       gff3 => 'http://sequenceontology.org/gff3/',
       gff3a => 'http://sequenceontology.org/gff3/attributes/',  # TEMP area
       oboInOwl => '"http://www.geneontology.org/formats/oboInOwl#',
       dc => 'http://purl.org/dc/terms/',
       part_of => 'http://purl.obolibrary.org/obo/BFO_0000050',   # TODO - check FALDO standard. RO vs SIO.
       sio => 'http://semanticscience.org/resource/SIO_',
       uniprotcore => 'http://purl.uniprot.org/core/',
      );

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

    my $strand =
        $f->strand ? $f->strand == 1 ? 'ForwardStrandPosition' : 'NegativeStrandPosition' 
        : 'BothStrandsPosition';
    $strand = "faldo:$strand";

    my $reference = $self->_id_to_uri($f->seq_id);

    # location object. In bioperl a feature can be split over multiple locations, or it
    # can have multiple locations. For now we assume a basic GFF like model where every
    # feature has a single location, so we create a single location object
    my @locs = ($self->_skolem($uri, 'loc'));
    foreach my $loc (@locs) {
        $self->_write_triple($uri, 'faldo:location', $loc);
        $self->_write_triple($loc, 'rdf:type', 'faldo:Region');

        # begin and end
        my $begin_obj = $self->_skolem($loc, 'b');
        my $end_obj = $self->_skolem($loc, 'e');
        $self->_write_triple($loc, 'faldo:begin', $begin_obj);
        $self->_write_triple($loc, 'faldo:end', $end_obj);
        $self->_write_triple($begin_obj, 'rdf:type', 'faldo:ExactPosition');
        $self->_write_triple($begin_obj, 'rdf:type', $strand);
        $self->_write_triple($begin_obj, 'faldo:position', _xsd( int=> $f->start ));
        $self->_write_triple($begin_obj, 'faldo:reference', $reference);
        $self->_write_triple($end_obj, 'rdf:type', 'faldo:ExactPosition');
        $self->_write_triple($end_obj, 'rdf:type', $strand);
        $self->_write_triple($end_obj, 'faldo:position', _xsd( int=> $f->end ));
        $self->_write_triple($end_obj, 'faldo:reference', $reference);
    }
        
    #source => $f->source_tag,

    if ($parent) {
        my ($parent_id) = $parent->get_tag_values('ID');
        if (!$parent_id) {
            $self->throw("feature ".($parent->display_name||'(no display name)')." has subfeatures but no 'ID' tag, cannot write as FALDO3");
        }
        $self->_write_triple($uri, 'part_of:', $self->_id_to_uri($parent_id));
    }

    # translate col9 (attributes) into triples
    foreach my $tag ( $f->get_all_tags ) {
        my ($p, $typefunc) = $self->_get_property_and_typefunc_for_tag($tag);

        # Note: order info will be lost; in future we may want a metadata tag for ordered properties
        my @vals = $f->get_tag_values( $tag );
        foreach my $v (@vals) {
            
            if ($tag eq 'Ontology_term') {
                $v =~ s/^(\S+):/obo:$1_/;
                # TODO - coordinate on model; avoid falling into OWL-Full
                $self->_write_triple($uri, 
                                     'uniprotcore:classifiedWith', 
                                     $v);
                next;
                                     
            }

            my $v_concrete = $v;
            if ($typefunc) {
                $v_concrete = $typefunc->($v);
            }
            else {
                if ($v =~ m@^[\-]?(\d+)@) {
                    $v_concrete = _xsd( int => $v );
                }
                elsif ($v =~ m@^[\-]?(\d+)\.\d+@) {
                    $v_concrete = _xsd( float => $v );
                }
                else {
                    $v_concrete = _quotify($v);
                }
            }
            $self->_write_triple($uri, $p, $v_concrete);
        }
    }


    $self->_print("\n");    
}


sub _get_property_and_typefunc_for_tag {
    my ($self, $tag) = @_;
    my $p = Bio::FeatureIO::OntologyMapper->get_tagmap->{lc($tag)};
    if ($p) {
        return ($p, \&_quotify );
    }
    return ("gff3a:$tag");
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

sub _quotify {
    my ($v) = @_;
    $v =~ s/\"/\\\"/g; # TODO - check this against TTL grammar
    return '"'.$v.'"';
}

sub _xsd {
    my ($t,$v) = @_;
    $v =~ s/\"/\\\"/g; # TODO - check this against TTL grammar
    return _quotify($v).'^^xsd:'.$t;
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

For more information in FALDO, see

  http://github.com/JervenBolleman/FALDO

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

 Chrus Mungall, <cjm at bioperl dot org>

Based on gff.pm, by:

 Chris Fields, <cjfields at bioperl dot org>
 Robert Buels, <rbuels at cpan dot org>
 Allen Day, <allenday@ucla.edu>


=head1 TODO

 - Fuzzy positions
 - Agree on URI naming conventions for injected nodes (location, begin, end)
 - Agree on vocabulary for gff3 attributes
 - Agree on model for Ontology_term

More broadly: try and standardize wider set of attributes used in GFF3, provide standard RDF vocab or OBO library mappings

=cut
    
