=pod

=head1 NAME

Bio::FeatureIO::bed - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

http://www.genome.ucsc.edu/goldenPath/help/customTrack.html#BED

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to
the Bioperl mailing list.  Your participation is much appreciated.

  bioperl-l@bioperl.org              - General discussion
  http://bioperl.org/MailList.shtml  - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
of the bugs and their resolution. Bug reports can be submitted via
the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

Describe contact details here

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::FeatureIO::bed;

use strict;
use base qw(Bio::FeatureIO);
use Bio::SeqFeature::Annotated;
use Bio::OntologyIO;

=head2 _initialize

 Title   : _initialize
 Function: initializes BED for reading/writing (currently write-only)
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
  my($self,%arg) = @_;

  $self->SUPER::_initialize(%arg);

  $self->name($arg{-name} || scalar(localtime()));
  $self->description($arg{-description} || scalar(localtime()));
  $self->use_score($arg{-use_score} || 0);

  $self->_print(sprintf('track name="%s" description="%s" useScore=%d',
                        $self->name,
                        $self->description,
                        $self->use_score ? 1 : 0
                       )
               );
}

=head2 use_score

 Title   : use_score
 Usage   : $obj->use_score($newval)
 Function: should score be used to adjust feature color when rendering?  set to true if so.
 Example : 
 Returns : value of use_score (a scalar)
 Args    : on set, new value (a scalar or undef, optional)


=cut

sub use_score{
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

sub name{
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

sub description{
    my $self = shift;

    return $self->{'description'} = shift if @_;
    return $self->{'description'};
}


sub write_feature {
  my($self,$feature) = @_;
  $self->throw("only Bio::SeqFeature::Annotated objects are writeable") unless $feature->isa('Bio::SeqFeature::Annotated');

  my $chrom       = $feature->seq_id || '';
  my $chrom_start = $feature->start  || 0;
  my $chrom_end   = $feature->stop   || 0;

  #try to make a reasonable name
  my $name        = undef;
  if(my @v = ($feature->annotation->get_Annotations('Name'))){
    $name = $v[0];
    $self->warn("only using first of feature's multiple names: ".join ',', map {$_->value} @v) if scalar(@v) > 1;
  } elsif(my @v = ($feature->annotation->get_Annotations('ID'))){
    $name = $v[0];
    $self->warn("only using first of feature's multiple IDs: ".join ',', map {$_->value} @v) if scalar(@v) > 1;
  } else {
    $name = 'anonymous';
  }

  my $score = $feature->score || 0;
  my $strand = $feature->strand == 0 ? '-' : '+'; #default to +
  my $thick_start = '';  #not implemented, used for CDS
  my $thick_end = '';    #not implemented, used for CDS
  my $reserved = 0;
  my $block_count = '';  #not implemented, used for sub features
  my $block_sizes = '';  #not implemented, used for sub features
  my $block_starts = ''; #not implemented, used for sub features

  $self->_print(join("\t",($chrom,$chrom_start,$chrom_end,$name,$score,$strand,$thick_start,$thick_end,$reserved,$block_count,$block_sizes, $block_starts)));
  $self->write_feature($_) foreach $feature->get_SeqFeatures();
}

sub next_feature {
  shift->throw_not_implemented();
}

1;
