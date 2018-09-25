package Bio::FeatureIO::ptt;

use utf8;
use strict;
use warnings;
use 5.010;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Generic;
use Data::Dumper;

# ABSTRACT: read/write features in PTT format
# AUTHOR:   Torsten Seemann <torsten.seemann@infotech.monash.edu.au>
# OWNER:    Torsten Seemann
# LICENSE:  Perl_5

# CONTRIBUTOR: based on bed.pm and gff.pm by Allen Day

=pod

=head1 SYNOPSIS

 # read features
 my $fin = Bio::FeatureIO->new(-file=>'genes.ptt', -format=>'ptt');
 my @cds;
 while (my $f = $fin->next_feature) {
   push @cds, $f if $f->strand > 0;
 }

 # write features (NOT IMPLEMENTED)
 my $fout = Bio::FeatureIO->new(-fh=>\*STDOUT, -format=>'ptt');
 for my $f (@cds) {
   $fout->write_feature($f);
 }

=head1 DESCRIPTION

The PTT file format is a table of protein features.
It is used mainly by NCBI who produce PTT files for
all their published genomes found in L<ftp://ftp.ncbi.nih.gov/genomes/>.
It has the following format:

=over 4

=item Line 1

Description of sequence to which the features belong
 eg. "Leptospira interrogans chromosome II, complete sequence - 0..358943"

It is usually equivalent to the DEFINITION line of a Genbank file,
with the length of the sequence appended. It is unclear why "0" is
used as a starting range, it should be "1".

=item Line 2

Number of feature lines in the table
 eg. "367 proteins"

=item Line 3

Column headers, tab separated
 eg. "Location  Strand  Length  PID Gene  Synonym Code  COG Product"

 Location : "begin..end" span of feature
 Strand   : "+" or "-"
 Length   : number of amino acids excluding the stop codon
 PID      : analogous to Genbank /db_xref="GI:xxxxxxxxx"
 Gene     : analogous to Genbank /gene="xxxx"
 Synonym  : analogous to Genbank /locus_tag="xxxx"
 Synonym  : analogous to Genbank /locus_tag="xxxx"
 COG      : CDD COG code with COG letter categories appended
 Product  : analogous to Genbank /product="xxxx"

=item Line 4 onwards

Feature lines, nine columns, tab separated, "-" used for empty fields
 eg. "2491..3423  + 310 24217063  metF  LB002 - COG0685E  5,10-methylenetetrahydrofolate reductase"


=back

=cut


# map tab-separated column number to field name
#our %NAME_OF = (
#    0 => 'Location',
#    1 => 'Strand',
#    2 => 'Length',
#    3 => 'PID',
#    4 => 'Gene',
#    5 => 'Synonym',
#    6 => 'Code',
#    7 => 'COG',
#    8 => 'Product',
#);
#our $NUM_COL = 9;

=head2 _initialize

 Title   : _initialize
 Function: Reading? parses the header of the input
           Writing?

=cut

#sub _initialize {
#    my ( $self, %arg ) = @_;
#
#    $self->SUPER::_initialize(%arg);
#
#    if ( $self->mode eq 'r' ) {
#
#        # Line 1
#        my $desc = $self->_readline();
#        chomp $desc;
#        $self->description($desc);
#
#        # Line 2
#        my $line = $self->_readline();
#        $line =~ m/^(\d+) proteins/ or $self->throw("Invalid protein count");
#        $self->protein_count($1);
#
#        # Line 3
#        $self->_readline();
#    }
#}

=head2 next_feature

 Title   : next_feature
 Usage   : $io->next_feature()
 Function: read the next feature from the PTT file
 Example :
 Args    :
 Returns : Bio::SeqFeatureI object

=cut

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
    my $dataset;
    while (defined($_ = $self->_readline)) {
        if (/^(\d+)..(\d+)\t/) {
            chomp;
            my (%feat, %tags);
            my @data = split("\t",$_);
            (@feat{qw(-location_string -strand -length -primary_id
                  -display_name)}, @tags{qw(Alias code cog product)}) = @data;
            for my $k (keys %tags) {
                delete $tags{$k} if $tags{$k} eq '-';
            }
            @feat{qw(-start -end)} = ($1, $2);
            $feat{-tag} = \%tags;
            @{$dataset}{qw(MODE DATA)} = ('feature', \%feat);
        } else {
            chomp;
            @{$dataset}{qw(MODE DATA)} = ('header', $_);
        }
        return $dataset;
    }
}

=head2 write_feature

 Title   : write_feature
 Usage   : $io->write_feature($feature)
 Function: write a Bio::SeqFeature object in PTT format with no header
 Example :
 Args    : Bio::SeqFeature object
 Returns :

=cut

sub write_feature {
  my($self,$feat) = @_;

  # Example, with header:
  # Location  Strand  Length  PID Gene  Synonym Code  COG Product
  # 190..255  + 21  16763391  thrL  STM0001 - - thr operon leader peptide

  $self->throw("Only Bio::SeqFeature::Generic or Bio::SeqFeature::Annotated objects are writeable")
  unless ( $feat->isa('Bio::SeqFeature::Generic') || $feat->isa('Bio::SeqFeature::Annotated') );

  # Default
  my ($len,$pid,$gene,$synonym,$code,$cog,$product) = qw(- - - - - - -);

  my $start = $feat->start;
  my $end   = $feat->end;
  my $loc   = "$start..$end";

  my $strand = $feat->strand == 1 ? '+' : '-';

  $len = int(($end - $start)/3);

  $product = join ' ',$feat->get_tag_values("product")     if ($feat->has_tag("product"));
  $pid     = join ' ',$feat->get_tag_values("protein_id")  if ($feat->has_tag("protein_id"));
  $code    = join ' ',$feat->get_tag_values("codon_start") if ($feat->has_tag("codon_start"));

  $self->_print(join("\t",($loc,$strand,$len,$pid,$gene,$synonym,$code,$cog,$product)) . "\n");

  $self->write_feature($_) foreach $feat->get_SeqFeatures();
}

=head2 description

 Title   : description
 Usage   : $obj->description($newval)
 Function: set/get the PTT file description for/from line one
 Example :
 Returns : value of description (a scalar)
 Args    : on set, new value (a scalar or undef, optional)

=cut

sub description {
    my $self = shift;
    $self->throw_not_implemented;
    # the above will just call the handler, which can handle persistent data
    # between calls
}

=head2 protein_count

 Title   : protein_count
 Usage   : $obj->protein_count($newval)
 Function: set/get the PTT protein count for/from line two
 Example :
 Args    : on set, new value (a scalar or undef, optional)
 Returns : value of protein_count (a scalar)

=cut

sub protein_count {
    my $self = shift;
    $self->throw_not_implemented;
    # the above will just call the handler, which can handle persistent data
    # between calls
}

1;
