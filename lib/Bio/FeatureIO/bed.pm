package Bio::FeatureIO::bed;

use utf8;
use strict;
use warnings;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Annotated;
use Bio::Annotation::SimpleValue;
use Scalar::Util qw(looks_like_number);

# ABSTRACT: read/write features from UCSC BED format
# AUTHOR:   Allen Day <allenday at ucla.edu>
# OWNER:    Allen Day
# LICENSE:  Perl_5

# CONTRIBUTOR: Sendu Bala <bix@sendu.me.uk>

=head2 _initialize

 Title   : _initialize
 Function: initializes BED for reading/writing
 Args    : all optional:
           name          description
           ----------------------------------------------------------
           -name         the name for the BED track, stored in header
                         name defaults to localtime()
           -description  the description for the BED track, stored in
                         header.  defaults to localtime().
           -use_score    whether or not the score attribute of
                         features should be used when rendering them.
                         the higher the score the darker the color.
                         defaults to 0 (false)



=cut

sub _initialize {
    my ( $self, %arg ) = @_;

    $self->SUPER::_initialize(%arg);

    $self->name( $arg{-name}               || scalar( localtime() ) );
    $self->description( $arg{-description} || scalar( localtime() ) );
    $self->use_score( $arg{-use_score}     || 0 );

    $self->_print(
        sprintf( 'track name="%s" description="%s" useScore=%d',
            $self->name, $self->description, $self->use_score ? 1 : 0 )
          . "\n"
    ) if $self->mode eq 'w';
}

=head2 use_score

 Title   : use_score
 Usage   : $obj->use_score($newval)
 Function: should score be used to adjust feature color when rendering?  set to true if so.
 Example :
 Returns : value of use_score (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub use_score {
    my $self = shift;

    return $self->{'use_score'} = shift if @_;
    return $self->{'use_score'};
}

=head2 name

 Title   : name
 Usage   : $obj->name($newval)
 Function: name of BED track
 Example :
 Returns : value of name (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub name {
    my $self = shift;

    return $self->{'name'} = shift if @_;
    return $self->{'name'};
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: description of BED track
 Example :
 Returns : value of description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub description {
    my $self = shift;

    return $self->{'description'} = shift if @_;
    return $self->{'description'};
}

sub next_feature {
    my $self = shift;
    while (my $ds = $self->next_dataset) {
        # leave it to the handler to decide when a feature is returned
        while (my $object = $self->handler->data_handler($ds)) {
            return $object if $object->isa('Bio::SeqFeatureI');
        }
    }
}

sub next_dataset {
    my $self = shift;
    my $mode = 'feature';
    while (my $line = $self->_readline) {
        my $dataset;
        next if $line =~ /^\s*$/;
        chomp $line;
        if ($line !~ /\t/) {
            $mode = 'track_definition';
            @{$dataset}{qw(MODE DATA)} = ($mode, {DATA => $line});
        } else {
            my (%feat, %tags);
            (@feat{qw(-seq_id -start -end -display_name -score -strand)},
             @tags{qw(thickstart thickend itemrgb blockct blocksizes blockstart)}) =
                    split(/\s+/, $line, 12);
            $feat{-strand} ||= '.'; # unknown is probably a good default
            $feat{-start} += 1;     # convert to 1-based coordinates
            $feat{-tag} = \%tags;
            @{$dataset}{qw(MODE DATA)} = ($mode, \%feat);
        }
        return $dataset;
    }
}

sub write_feature {
    my ( $self, $feature ) = @_;
    $self->throw_not_implemented;
    $self->throw("only Bio::SeqFeature::Annotated objects are writeable")
      unless $feature->isa('Bio::SeqFeature::Annotated');

    my $chrom = $feature->seq_id || '';
    my $chrom_start = $feature->start
      || 0;    # output start is supposed to be 0-based
    my $chrom_end = ( $feature->end + 1 )
      || 1;    # output end is supposed to not be part of the feature

    #try to make a reasonable name
    my $name = undef;
    my @v;
    if ( @v = ( $feature->annotation->get_Annotations('Name') ) ) {
        $name = $v[0];
        $self->warn(
            "only using first of feature's multiple names: " . join ',',
            map { $_->value } @v )
          if scalar(@v) > 1;
    }
    elsif ( @v = ( $feature->annotation->get_Annotations('ID') ) ) {
        $name = $v[0];
        $self->warn( "only using first of feature's multiple IDs: " . join ',',
            map { $_->value } @v )
          if scalar(@v) > 1;
    }
    else {
        $name = 'anonymous';
    }

    if ( ref($name) ) {
        $name = $name->value;
    }
    if ( ref($chrom) ) {
        $chrom = $chrom->value;
    }

    my $score = $feature->score || 0;
    my $strand = $feature->strand == 0 ? '-' : '+';    #default to +
    my $thick_start  = '';    #not implemented, used for CDS
    my $thick_end    = '';    #not implemented, used for CDS
    my $reserved     = 0;
    my $block_count  = '';    #not implemented, used for sub features
    my $block_sizes  = '';    #not implemented, used for sub features
    my $block_starts = '';    #not implemented, used for sub features

    $self->_print(
        join(
            "\t",
            (
                $chrom,    $chrom_start, $chrom_end,   $name,
                $score,    $strand,      $thick_start, $thick_end,
                $reserved, $block_count, $block_sizes, $block_starts
            )
          )
          . "\n"
    );
    $self->write_feature($_) foreach $feature->get_SeqFeatures();
}

1;

__END__

=pod

=head1 SYNOPSIS

  my $in = Bio::FeatureIO->new(-format => 'bed', -file => 'file.bed');
  for my $feat ($in->next_feature) {
    # do something with $feat (a Bio::SeqFeature::Annotated object)
  }

  my $out = Bio::FeatureIO->new(-format=>'bed');
  for my $feat ($seq->get_seqFeatures) {
    $out->write_feature($feat);
  }

=head1 DESCRIPTION

See L<http://www.genome.ucsc.edu/goldenPath/help/customTrack.html#BED>.

Currently for read and write only the first 6 fields (chr, start, end, name,
score, strand) are supported.

=cut
