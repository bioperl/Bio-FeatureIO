package Bio::FeatureIO::faldo;

use utf8;
use strict;
use warnings;
use 5.010;
use base qw(Bio::FeatureIO);
use Bio::FeatureIO::OntologyMapper;

# ABSTRACT: read/write FALDO feature files
# AUTHOR:   Chrus Mungall <cjm at bioperl dot org>
# OWNER:    Chrus Mungall
# LICENSE:  Perl_5

##Based on gff.pm, by:
##
## Chris Fields, <cjfields at bioperl dot org>
## Robert Buels, <rbuels at cpan dot org>
## Allen Day, <allenday@ucla.edu>

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
        $self->_write_triple($begin_obj, 'faldo:position', _xsd( int=> $f->strand >= 0 ? $f->start : $f->end ));
        $self->_write_triple($begin_obj, 'faldo:reference', $reference);
        $self->_write_triple($end_obj, 'rdf:type', 'faldo:ExactPosition');
        $self->_write_triple($end_obj, 'rdf:type', $strand);
        $self->_write_triple($end_obj, 'faldo:position', _xsd( int=> $f->strand >= 0 ?  $f->end : $f->start ));
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

sub _so {
    my ($self, $t, $is_full_uri) = @_;

    unless ($is_full_uri) {
        if ($self->{_prefixmap_written}->{$t}) {
            return "$t:";
        }
    }
    my $so_id = $self->{_typemap}->{$t};
    if ($so_id) {
        if ($is_full_uri) {
            $so_id =~ s@^SO:@http://purl.obolibrary.org/obo/SO_@;
        }
        else {
            $so_id =~ s/^SO:/obo:SO_/;
        }

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
        $self->_write_prefix_declaration($t, $self->_so($t, 1));
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

=head1 SYNOPSIS

  my $feature; #get a Bio::SeqFeatureI somehow
  my $featureOut = Bio::FeatureIO->new(
    -format => 'faldo',
    -fh => \*STDOUT,
  );
  $featureOut->write_feature($feature);

=head1 DESCRIPTION

This writer produces turtle-formatted RDF following the FALDO schema
from bioperl objects.

For more information in FALDO, see

  http://github.com/JervenBolleman/FALDO

Note: Parsing is not yet implemented. It's not clear if bioperl is the
appropriate vehicle for this.

=head2 IMPLEMENTATION

=head3 Turtle export

Rather than introduce a dependency on an RDF library, we export turtle
directly. This give us some control over the output and allows us to
make the output more human-readable.

For example, the header declares a prefix for every SO type:

  @prefix exon: <http://purl.obolibrary.org/obo/SO_0000147> .
  @prefix intron: <http://purl.obolibrary.org/obo/SO_0000188> .
  @prefix mRNA: <http://purl.obolibrary.org/obo/SO_0000234> .

This allows us to use the human-readable prefix in place of the URI,
for example:

  <http://bioperl.org/features/mRNA:Solyc03g123530.2.1> rdf:type mRNA: .

=head3 Generation of URIs

Rather than using RDF blank nodes, we auto-generate/skolemize URIs for
all entities.

 - Feature URIs are generated by prefixing the feature ID with a URI prefix (default bioperl.org)
 - Location URIs are generated by suffixing __loc onto the feature URI (assuming one location per feature)
 - Begin/end URIs are generated by suffixing __b or __e onto the Location URI

=head1 TODO

 - Fuzzy positions
 - Agree on URI naming conventions for injected nodes (location, begin, end)
 - Agree on vocabulary for gff3 attributes
 - Agree on model for Ontology_term

More broadly: try and standardize wider set of attributes used in GFF3, provide standard RDF vocab or OBO library mappings

=cut
