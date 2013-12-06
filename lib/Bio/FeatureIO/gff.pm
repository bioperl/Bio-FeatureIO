package Bio::FeatureIO::gff;

use strict;
use warnings;
use 5.010;
use base qw(Bio::FeatureIO);
# TODO: work on integrating more of Rob's work here. Fewer points of failure,
# the better...

use Bio::GFF3::LowLevel qw(
    gff3_format_feature
    gff3_escape
    gff3_unescape
);

sub _initialize {
    my ($self, @args) = @_;
    $self->SUPER::_initialize(@args);

    my %args = @args;
    $self->ignore_seq_region($args{-ignore_seq_region} || 0);

    my ($version) = $self->_rearrange([qw(VERSION)], @args);
    $version ||= 3;
    $self->version($version);
    # set globals using config if present, then defaults based on version
}

# raw feature stream; returned features are as-is, may be modified post-return
sub next_feature {
    my $self = shift;
    return if $self->fasta_mode;
    DATASET:
    while (my $ds = $self->next_dataset) {
        # leave it to the handler to decide when a feature is returned
        while (my $object = $self->handler->data_handler($ds)) {
            return $object if $object->isa('Bio::SeqFeatureI');
            # when a SeqIO is returned, the features are done
            if ($object->isa('Bio::SeqIO')) {
                $self->seqio($object);
                return;
            }
        }
    }
}

sub next_dataset {
    my $self = shift;

    # Note: this simply allows for versions of escape/unescape functions TODO:
    # needs to be made more stable and specifically versioned (only GFF3 seems
    # to be supported)
    
    state $unescape_func = {
        3 => \&gff3_unescape,
      }->{$self->version}
      || sub { $_[0] };

    state $escape_func = {
        3 => \&gff3_escape,
      }->{$self->version}
      || sub { $_[0] };

    # In Windows, text files have '\r\n' as line separator, but when reading in
    # text mode Perl will only show the '\n'. This means that for a line "ABC\r\n",
    # "length $_" will report 4 although the line is 5 bytes in length.
    # We assume that all lines have the same line separator and only read current line.
    my $fh         = $self->_fh;
    my $init_pos   = tell($fh);
    my $init_line  = $.;
    my $curr_line  = <$fh>;
    my $pos_diff   = tell($fh) - $init_pos;
    # If position difference is 0, cursor is already at the end of the file
    my $correction = ($pos_diff > 0) ? ($pos_diff - length $curr_line) : 0;
    $fh->input_line_number($init_line); # Rewind line number $.
    seek $fh, $init_pos, 0;             # Rewind position to proceed to read the file

    local $/ = "\n";
    my $dataset;
    my $len = 0;
    GFFLINE:
    while (my $line = $self->_readline) {
        $len += CORE::length($line) + $correction;
        for ($line) {
            if (/^\s*$/) {  next GFFLINE  } # blank lines
            elsif (/^(\#{1,2})\s*(\S+)\s*([^\n]+)?$/) { # comments and directives
                if (length($1) == 1) {
                    chomp $line;
                    @{$dataset}{qw(MODE DATA)} = ('comment', {DATA => $line});
                } else {
                    $self->{mode} = 'directive';
                    @{$dataset}{qw(MODE DATA)} =
                        ('directive', $self->directive($2, $3));
                }
            }
            elsif (/^>/) {          # sequence
                chomp $line;
                @{$dataset}{qw(MODE DATA)} =
                    ('sequence', {'sequence-header' =>  $line});
                $self->{mode} = 'sequence';
            }
            elsif (/(?:\t[^\t]+){8}/)  {
                chomp $line;
                $self->{mode} = $dataset->{MODE} = 'feature';
                my (%feat, %tags, $attstr);
                # validate here?
                (@feat{qw(-seq_id -source -primary_tag -start -end
                       -score -strand -phase)}, $attstr) =
                    map {$_ eq '.' ? undef : $_ } split( "\t" ,$line);

                for my $kv (split(/\s*;\s*/, $attstr)) {
                    my ($key, $rest) = split( '=', $kv, 2);
                    my @vals = $rest ? map { $unescape_func->($_) } split(',',$rest)
                        : ();
                    push @{$tags{$key} ||= [] }, @vals;
                }
                $feat{-tag} = \%tags;
                $dataset->{DATA} = \%feat;
            }
            else {
                if ($self->{mode} eq 'sequence') {
                    chomp $line;
                    @{$dataset}{qw(MODE DATA)} =
                        ('sequence', {sequence => $line});
                } else {
                    # anything else should be sequence, but there should be some
                    # kind of directive to change the mode or a typical FASTA
                    # header should be found; if not, die
                    $self->throw("Unknown line: $line, parser was in mode ".
                                 $self->{mode});
                }
            }
        }
        if ($dataset) {
            @$dataset{qw(START LENGTH)} = ($self->{stream_start}, $len);
            $self->{stream_start} += $len;
            return $dataset;
        }
        return;
    }
}

sub directive {
    my ($self, $directive, $rest) = @_;
    $rest ||= '';
    my %data;

    # note this is allowing for expansion of additional directives
    for ($directive) {
        # validate here?
        if ($_ eq 'sequence-region') {
            next if $self->ignore_seq_region();
            @data{qw(type id start end)} =
                ('sequence-region', split(/\s+/, $rest));
        }
        elsif ($_ eq 'genome-build') {
            @data{qw(type source buildname)} = ($directive, split(/\s+/, $rest));
        }
        elsif ($_ eq '#') {
            $data{type} = 'resolve-references';
        }
        elsif ($_ eq 'FASTA') {
            $data{type} = 'sequence';
        }
        else {
            @data{qw(type data)} = ($directive, $rest);
        }
    }
    \%data;
}

=head1 ignore_seq_region

Set this flag to keep FeatureIO from returning
a feature for a ##sequence-region directive

=cut

sub ignore_seq_region {
  my($self,$val) = @_;
  $self->{'ignore_seq_region'} = $val if defined($val);
  return $self->{'ignore_seq_region'};
}

# TODO: this gets into the handler internals a bit too much, and I think the
# problem is orthogonal to the feature stream (in other words, maybe we simply
# pass features along and plugin(s) will help resolve feature groups)

sub next_feature_group {
    my $self = shift;
    return if $self->fasta_mode; # end of features in the file
    my %seen_ids;
    my @all_feats;
    my @toplevel_feats;

    while (my $ds = $self->next_dataset) {
        my $object = $self->handler->data_handler($ds);
        if ($object && $object->isa('Bio::SeqIO')) {
            $self->seqio($object);
            last;
        }

        if( $self->handler->resolve_references ) {
            if( @toplevel_feats ) {
                last;
            } else {
                next;
            }
        }

        next unless $object;

        if ($object->has_tag('ID')) {
            my ($id) = $object->get_tag_values('ID');
            $self->throw("Oops! ID $id exists more than once in your file!")
                if (exists($seen_ids{$id}));
            $seen_ids{$id} = $object;
            push @all_feats, $object;
            push @toplevel_feats, $object if !$object->has_tag('Parent');
        }
        if ($object->has_tag('Parent')) {
            my @parents = $object->get_tag_values('Parent');
            for my $parent_id (@parents) {
                if (exists $seen_ids{$parent_id}) {
                    $seen_ids{$parent_id}->add_SeqFeature($object);
                } else {
                    $self->throw("Parent with ID $parent_id not found!");
                }
            }
        }
    }
    return @toplevel_feats;
}

sub next_seq() {
    my $self = shift;
    return undef unless $self->fasta_mode();
    return $self->seqio->next_seq();
}

=head2 write_feature()

 Usage   : $featureio->write_feature( Bio::SeqFeature::Annotated->new(...) );
 Function: writes a feature in GFF format.  the GFF version used is
           governed by the '-version' argument passed to Bio::FeatureIO->new(),
           and defaults to GFF version 3.
 Returns : nothing meaningful
 Args    : a Bio::SeqFeature::Annotated object.

=cut

sub write_feature {
    my($self,$feature) = @_;

    unless( $self->{wrote_gff_version} ) {
        $self->_print("##gff-version ".$self->version."\n");
        $self->{wrote_gff_version} = 1;
    }

    unless( $feature ) {
        $self->throw("gff.pm cannot write_feature unless you give a feature to write.\n");
    }

    my $funcname = '_write_feature_'.$self->version;
    $self->can( $funcname )
       or $self->throw( 'writing not implemented for GFF version '.$self->version );
    { no strict 'refs';
      return $self->$funcname( $feature );
    }
}

sub _write_feature_3 {
    my ( $self, $feature, $parent_feature ) = @_;

    my $str = gff3_format_feature( $self->_gff3_lowlevel_hashref( $feature, $parent_feature ))
        # TODO: add some more info about the feature to this error message
        or $self->throw( 'failed to format feature for writing' );
    $self->_print( $str );

    $self->_write_feature_3( $_, $feature ) for $feature->get_SeqFeatures;

    $self->_print( "###\n" ) unless $parent_feature;
}

sub _gff3_lowlevel_hashref {
    my ( $self, $f, $parent ) = @_;

    my @tags = $f->get_all_tags;
    if( $f->can('phase') ) {
        @tags = grep $_ ne 'phase', @tags;
    }

    if( $parent && ! $parent->has_tag('ID') ) {
        $self->throw("feature ".($parent->display_name||'(no display name)')." has subfeatures but no 'ID' tag, cannot write as GFF3");
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

################################################################################

=head1 ACCESSORS

=cut

=head2 fasta_mode()

 Usage   : $obj->fasta_mode($newval)
 Function: Indicates whether or not the parser is currently processing the FASTA
           section of the GFF stream
 Example : 
 Returns : Value of fasta_mode (a scalar)
 Args    : None

=cut

sub fasta_mode {
    return shift->handler->fasta_mode;
}

=head2 seqio()

 Usage   : $obj->seqio($newval)
 Function: get/set a Bio::SeqIO instance to handle the GFF3 ##FASTA section.
 Returns : a Bio::SeqIO object or undef
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub seqio {
  my($self,$val) = @_;
  if (defined $val) {
    $self->{'seqio'} = $val;
  } else {
    # Cannot get seqio before we've reached the ##FASTA section
    return undef unless $self->fasta_mode();
    if (not defined $self->{'seqio'}) {
      # Initialize Bio::SeqIO instance
      $self->{'seqio'} = Bio::SeqIO->new(-format => 'fasta', -fh => $self->_fh());
    }
  }
  return $self->{'seqio'};
}

=head2 sequence_region()

 Usage   :
 Function: ###FIXME
 Returns : 
 Args    :

=cut

sub sequence_region {
    shift->throw_not_implemented;
#    my ($self,$k,$v) = @_;
#    if(defined($k) && defined($v)){
#        $self->{'sequence_region'}{$k} = $v;
#        return $v;
#    }
#    elsif(defined($k)){
#        return $self->{'sequence_region'}{$k};
#    }
#    else {
#        return;
#    }
}


=head2 so()

 Usage   : $obj->so($newval)
 Function: holds a Sequence Ontology instance
 Returns : value of so (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub so {
    shift->throw_not_implemented;
    #my $self = shift;
    #my $val = shift;
    ####FIXME validate $val object's type
    #$self->{so} = $val if defined($val);
    #return $self->{so};
}

=head2 validate()

 Usage   : $obj->validate($newval)
 Function: true if encountered ontology terms in next_feature()
           mode should be validated.
 Returns : value of validate (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub validate {
    shift->throw_not_implemented;
    #my($self,$val) = @_;
    #$self->{'validate'} = $val if defined($val);
    #return $self->{'validate'};
}

=head2 version()

 Usage   : $obj->version($newval)
 Function: version of GFF to read/write.  valid values are 1, 2, 2.5, and 3.
 Returns : value of version (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

{
    my %VALID_VERSION = map {$_ => 1} (1, 2, 2.5, 3);

sub version {
    my $self = shift;
    my $val = shift;
    if(defined $val){
        $self->throw('Invalid GFF version.  Valid versions: '.join(' ', sort keys %VALID_VERSION))
            if !exists($VALID_VERSION{$val});
        $self->{'version'} = $val;
        #$self->set_config();
    }
    return $self->{'version'};
}

}

1;

__END__

=head1 NAME

Bio::FeatureIO::gff - read/write GFF feature files

=head1 SYNOPSIS

  my $feature; #get a Bio::SeqFeatureI somehow
  my $featureOut = Bio::FeatureIO->new(
    -format => 'gff',
    -version => 3,
    -fh => \*STDOUT,
    -validate_terms => 1, #boolean. validate ontology terms online?  default 0 (false).
  );
  $featureOut->write_feature($feature);

=head1 DESCRIPTION

This is a general purpose GFF parser.  The default is GFF version 3.

Data is passed as hash-refs, similar to a SAX-based data stream, but containing
chunks of related information.  A version of this is implemented in Bio::SeqIO
plugins gbdriver, embldriver, and swissdriver.

The key issue is defining specifically how bits are bundled and passed along to
the data handler. the other key point is that the start and length of the
specific chunk of data passed in is also passed along, primarily if one wanted
to create lazy feature collections.

The structure of the passed hash references is possibly in flux and shouldn't be
directly relied on; I plan on standardizing this for consistency. One
possibility is to standardize on BioPerl constructor attributes where possible,
which (in general) tends to mirror many natural data types such as features.

Maybe something like:

 # GFF directives
 $VAR = {
    'MODE'      => 'directive', # top level type
    'DATA'      => {
        'PRIMARY_TYPE'   => 'VERSION',  # second level (possible subtype)
        # what follows are specific to type/subtype pairings
        'VERSION'       => 3
    }
 };
 
 # features
 
 $VAR = {
    'TYPE'      => 'FEATURE',
    'DATA'      => {
        'PRIMARY_TYPE'   => 'SEQUENCE_FEATURE',
        'TYPE'          => 'gene',
        'SOURCE'        => 'GenBank',
        'START'         => 1,
        'END'           => 10000,
        'STRAND'        => -1,
        'PHASE'         => '.', # can also be left out 
        'SCORE'         => '.',
        'ATTRIBUTES'    => {
            'name'      => ['BRCA1'],
            'dbxref'    => [...],
        }
    }
 };
 
 Currently implemented:

 version         read?   write?
 ------------------------------
 GFF 1             N       N
 GFF 2             N       N
 GFF 2.5 (GTF)     N       Y
 GFF 3             Y       Y

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
 Function: reads a feature record from a GFF stream and returns it as an object.
 Returns : a Bio::SeqFeature::Annotated object
 Args    : N/A

=head2 next_feature_group

 Usage   : @feature_group = $stream->next_feature_group
 Function: Reads the next feature_group from $stream and returns it.

           Feature groups in GFF3 files are separated by '###' directives. The
           features in a group might form a hierarchical structure. The
           complete hierarchy of features is returned, i.e. the returned array
           represents only the top-level features.  Lower-level features can
           be accessed using the 'get_SeqFeatures' method recursively.

 Example : # getting the complete hierarchy of features in a GFF3 file
           my @toplevel_features;
           while (my @fg = $stream->next_feature_group) {
               push(@toplevel_features, @fg);
           }
 Returns : an array of Bio::SeqFeature::Annotated objects
 Args    : none

=head2 next_seq()

  Usage   : $featureio->next_seq( );
  Function: access the FASTA section (if any) at the end of the GFF stream. Note
            that this method will return undef before all the features in the GFF
            stream have been handled.
  Returns : a Bio::SeqI object or undef
  Args    : none

=cut
