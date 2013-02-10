# $Id$
#
# BioPerl module for Bio::FeatureIO
#
# Cared for by Allen Day <allenday@ucla.edu>
#
# Copyright Allen Day
#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code

=head1 NAME

Bio::FeatureIO - Handler for FeatureIO

=head1 SYNOPSIS

  use Bio::FeatureIO;

  #read from a file
  $in  = Bio::FeatureIO->new(-file => "my.gff" , -format => 'GFF');

  #read from a filehandle
  $in  = Bio::FeatureIO->new(-fh => \*GFF , -format => 'GFF');

  #read features already attached to a sequence
  my $seq = Bio::Seq->new(-seq => $seq , -format => 'features');

  #read new features for existing sequence
  my $seq = Bio::Seq->new(-seq => $seq , -format => 'Das');

  #write out features
  $out = Bio::FeatureIO->new(-file    => ">outputfilename" ,
                             -format  => 'GFF3' ,
                             -version => 3);

  while ( my $feature = $in->next_feature() ) {
    $out->write_feature($feature);
  }

=head1 DESCRIPTION

FIXME

=head1 SUPPORTED FORMATS

 name                      module
 --------------------------------
 BED                       bed.pm
 GFF                       gff.pm
 GTF                       gtf.pm

#Bio::SeqIO is a handler module for the formats in the SeqIO set (eg,
#Bio::SeqIO::fasta). It is the officially sanctioned way of getting at
#the format objects, which most people should use.
#
#The Bio::SeqIO system can be thought of like biological file handles.
#They are attached to filehandles with smart formatting rules (eg,
#genbank format, or EMBL format, or binary trace file format) and 
#can either read or write sequence objects (Bio::Seq objects, or
#more correctly, Bio::SeqI implementing objects, of which Bio::Seq is
#one such object). If you want to know what to do with a Bio::Seq
#object, read L<Bio::Seq>.
#
#The idea is that you request a stream object for a particular format.
#All the stream objects have a notion of an internal file that is read
#from or written to. A particular SeqIO object instance is configured
#for either input or output. A specific example of a stream object is
#the Bio::SeqIO::fasta object.
#
#Each stream object has functions
#
#   $stream->next_feature();
#
#and
#
#   $stream->write_feature($feature);
#

=head1 CONSTRUCTORS

=head2 Bio::FeatureIO-E<gt>new()

   $featureIO = Bio::FeatureIO->new(-file => 'filename',   -format=>$format);
   $featureIO = Bio::FeatureIO->new(-fh   => \*FILEHANDLE, -format=>$format);
   $featureIO = Bio::FeatureIO->new(-seq  => $seq,         -format=>$format);

The new() class method constructs a new Bio::FeatureIO object.  The
returned object can be used to retrieve or print Seq objects. new()
accepts the following parameters:

=over 4

=item -file

A file path to be opened for reading or writing.  The usual Perl
conventions apply:

   'file'       # open file for reading
   '>file'      # open file for writing
   '>>file'     # open file for appending
   '+<file'     # open file read/write
   'command |'  # open a pipe from the command
   '| command'  # open a pipe to the command

=item -fh

You may provide new() with a previously-opened filehandle.  For
example, to read from STDIN:

   $seqIO = Bio::FeatureIO->new(-fh => \*STDIN);

Note that you must pass filehandles as references to globs.

If neither a filehandle nor a filename is specified, then the module
will read from the @ARGV array or STDIN, using the familiar E<lt>E<gt>
semantics.

A string filehandle is handy if you want to modify the output in the
memory, before printing it out. The following program reads in EMBL
formatted entries from a file and prints them out in fasta format with
some HTML tags:

  use Bio::SeqIO;
  use IO::String;
  my $in  = Bio::SeqIO->new('-file' => "emblfile" , 
  			    '-format' => 'EMBL');
  while ( my $seq = $in->next_seq() ) {
      # the output handle is reset for every file
      my $stringio = IO::String->new($string);
      my $out = Bio::SeqIO->new('-fh' => $stringio,
  			        '-format' => 'fasta');
      # output goes into $string
      $out->write_seq($seq);
      # modify $string
      $string =~ s|(>)(\w+)|$1<font color="Red">$2</font>|g;
      # print into STDOUT
      print $string;
  }

=item -format

Specify the format of the file.  Supported formats include:


=item -flush

By default, all files (or filehandles) opened for writing sequences
will be flushed after each write_seq() (making the file immediately
usable).  If you don't need this facility and would like to marginally
improve the efficiency of writing multiple sequences to the same file
(or filehandle), pass the -flush option '0' or any other value that
evaluates as defined but false:

  my $f1 = new Bio::FeatureIO -file   => "<a.f1",
                              -format => "f1";
  my $f2 = new Bio::FeatureIO -file   => ">a.f2",
                              -format => "f2",
                              -flush  => 0; # go as fast as we can!

  while($feature = $f1->next_seq) { $f2->write_seq($feature) }

=back

=head2 Bio::FeatureIO-E<gt>newFh()

   $fh = Bio::FeatureIO->newFh(-fh   => \*FILEHANDLE, -format=>$format);
   $fh = Bio::FeatureIO->newFh(-format => $format);
   # etc.

This constructor behaves like new(), but returns a tied filehandle
rather than a Bio::FeatureIO object.  You can read sequences from this
object using the familiar E<lt>E<gt> operator, and write to it using
print().  The usual array and $_ semantics work.  For example, you can
read all sequence objects into an array like this:

  @sequences = <$fh>;

Other operations, such as read(), sysread(), write(), close(), and printf()
are not supported.

=head1 OBJECT METHODS

See below for more detailed summaries.  The main methods are:

=head2 $feature = $featureIO-E<gt>next_feature()

Fetch the next feature from the stream.

=head2 $featureIO-E<gt>write_feature($feature [,$another_feature,...])

Write the specified feature(s) to the stream.

=head2 TIEHANDLE(), READLINE(), PRINT()

These provide the tie interface.  See L<perltie> for more details.

=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this
and other Bioperl modules. Send your comments and suggestions preferably
to one of the Bioperl mailing lists.

Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/MailList.shtml      - About the mailing lists

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
 the bugs and their resolution.
 Bug reports can be submitted via the web:

  http://bugzilla.bioperl.org/

=head1 AUTHOR - Allen Day

Email allenday@ucla.edu

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

#' Let the code begin...

package Bio::FeatureIO;

use strict;
use vars qw(@ISA);

use Bio::Root::Root;
use Bio::Root::IO;
use Symbol();

@ISA = qw(Bio::Root::Root Bio::Root::IO);

=head2 new

 Title   : new
 Usage   : $stream = Bio::FeatureIO->new(-file => $filename, -format => 'Format')
 Function: Returns a new seqstream
 Returns : A Bio::FeatureIO stream initialised with the appropriate format
 Args    : Named parameters:
             -file => $filename
             -fh => filehandle to attach to
             -format => format

=cut

my $entry = 0;

sub new {
  my ($caller,@args) = @_;
  my $class = ref($caller) || $caller;

  # or do we want to call SUPER on an object if $caller is an
  # object?
  if( $class =~ /Bio::FeatureIO::(\S+)/ ) {

    my ($self) = $class->SUPER::new(@args);	
    $self->_initialize(@args);
    return $self;

  } else {

	my %param = @args;

	@param{ map { lc $_ } keys %param } = values %param; # lowercase keys
	my $format = $param{'-format'} ||
      $class->_guess_format( $param{-file} || $ARGV[0] );
	
	if( ! $format ) {
      if ($param{-file}) {
		$format = Bio::Tools::GuessSeqFormat->new(-file => $param{-file}||$ARGV[0] )->guess;
      } elsif ($param{-fh}) {
		$format = Bio::Tools::GuessSeqFormat->new(-fh => $param{-fh}||$ARGV[0] )->guess;
      }
	}
	$format = "\L$format";	# normalize capitalization to lower case
	return undef unless( $class->_load_format_module($format) );
	return "Bio::FeatureIO::$format"->new(@args);

  }
}

=head2 newFh

 Title   : newFh
 Usage   : $fh = Bio::FeatureIO->newFh(-file=>$filename,-format=>'Format')
 Function: does a new() followed by an fh()
 Example : $fh = Bio::FeatureIO->newFh(-file=>$filename,-format=>'Format')
           $feature = <$fh>;   # read a feature object
           print $fh $sequence; # write a feature object
 Returns : filehandle tied to the Bio::FeatureIO::Fh class
 Args    :

See L<Bio::FeatureIO::Fh>

=cut

sub newFh {
  my $class = shift;
  return unless my $self = $class->new(@_);
  return $self->fh;
}

=head2 fh

 Title   : fh
 Usage   : $obj->fh
 Function:
 Example : $fh = $obj->fh;      # make a tied filehandle
           $sequence = <$fh>;   # read a sequence object
           print $fh $sequence; # write a sequence object
 Returns : filehandle tied to Bio::FeatureIO class
 Args    : none

=cut


sub fh {
  my $self = shift;
  my $class = ref($self) || $self;
  my $s = Symbol::gensym;
  tie $$s,$class,$self;
  return $s;
}

# _initialize is chained for all FeatureIO classes

sub _initialize {
    my($self, %arg) = @_;

    # flush is initialized by the Root::IO init

    # initialize the IO part
    $self->seq($arg{-seq});
    $self->_initialize_io(%arg);
}

=head2 next_feature

 Title   : next_feature
 Usage   : $feature = stream->next_feature
 Function: Reads the next feature object from the stream and returns it.

           Certain driver modules may encounter entries in the stream
           that are either misformatted or that use syntax not yet
           understood by the driver. If such an incident is
           recoverable, e.g., by dismissing a feature of a feature
           table or some other non-mandatory part of an entry, the
           driver will issue a warning. In the case of a
           non-recoverable situation an exception will be thrown.  Do
           not assume that you can resume parsing the same stream
           after catching the exception. Note that you can always turn
           recoverable errors into exceptions by calling
           $stream->verbose(2).

 Returns : a Bio::SeqFeatureI feature object
 Args    : none

See L<Bio::Root::RootI>, L<Bio::SeqFeatureI>

=cut

sub next_feature {
   my ($self, $seq) = @_;
   $self->throw("Sorry, you cannot read from a generic Bio::FeatureIO object.");
}

=head2 write_feature

 Title   : write_feature
 Usage   : $stream->write_feature($feature)
 Function: writes the $feature object into the stream
 Returns : 1 for success and 0 for error
 Args    : Bio::SeqFeature object

=cut

sub write_feature {
    my ($self, $seq) = @_;
    if(ref($self) eq __PACKAGE__){
      $self->throw("Sorry, you cannot write to a generic Bio::FeatureIO object.");
    } else {
      $self->throw_not_implemented();
    }
}

=head2 _load_format_module

 Title   : _load_format_module
 Usage   : *INTERNAL FeatureIO stuff*
 Function: Loads up (like use) a module at run time on demand
 Example :
 Returns :
 Args    :

=cut

sub _load_format_module {
    my ($self, $format) = @_;
    my $class = ref($self) || $self;
    my $module = $class."::$format";#"Bio::Feature::" . $format;
    my $ok;

    eval {
	$ok = $self->_load_module($module);
    };
    if ( $@ ) {
    print STDERR <<END;
$self: $format cannot be found
Exception $@
For more information about the FeatureIO system please see the FeatureIO docs.
This includes ways of checking for formats at compile time, not run time
END
  ;
  }
  return $ok;
}

=head2 seq()

=over

=item Usage

  $obj->seq();        #get existing value

  $obj->seq($newval); #set new value

=item Function


=item Returns

value of seq (a Bio::SeqI)

=item Arguments

(optional) on set, a scalar

=back

=cut

sub seq {
  my $self = shift;
  my $val = shift;

  $self->{'seq'} = $val if defined($val);
  return $self->{'seq'};
}


=head2 _filehandle

 Title   : _filehandle
 Usage   : $obj->_filehandle($newval)
 Function: This method is deprecated. Call _fh() instead.
 Example :
 Returns : value of _filehandle
 Args    : newvalue (optional)


=cut

sub _filehandle {
    my ($self,@args) = @_;
    return $self->_fh(@args);
}

=head2 _guess_format

 Title   : _guess_format
 Usage   : $obj->_guess_format($filename)
 Function: guess format based on file suffix
 Example :
 Returns : guessed format of filename (lower case)
 Args    :
 Notes   : formats that _filehandle() will guess include FIXME...

=cut

sub _guess_format {
   my $class = shift;
   return unless $_ = shift;
   return 'fasta'   if /\.(fasta|fast|fas|seq|fa|fsa|nt|aa)$/i;
   return 'genbank' if /\.(gb|gbank|genbank|gbk|gbs)$/i;
   return 'scf'     if /\.scf$/i;
   return 'abi'     if /\.ab[i1]$/i;
   return 'alf'     if /\.alf$/i;
   return 'ctf'     if /\.ctf$/i;
   return 'ztr'     if /\.ztr$/i;
   return 'pln'     if /\.pln$/i;
   return 'exp'     if /\.exp$/i;
   return 'pir'     if /\.pir$/i;
   return 'embl'    if /\.(embl|ebl|emb|dat)$/i;
   return 'raw'     if /\.(txt)$/i;
   return 'gcg'     if /\.gcg$/i;
   return 'ace'     if /\.ace$/i;
   return 'bsml'    if /\.(bsm|bsml)$/i;
   return 'swiss'   if /\.(swiss|sp)$/i;
   return 'phd'     if /\.(phd|phred)$/i;
   return 'fastq'   if /\.fastq$/i;
}

sub DESTROY {
    my $self = shift;

    $self->close();
}

sub TIEHANDLE {
    my ($class,$val) = @_;
    return bless {'seqio' => $val}, $class;
}

sub READLINE {
  my $self = shift;
  return $self->{'seqio'}->next_seq() unless wantarray;
  my (@list, $obj);
  push @list, $obj while $obj = $self->{'seqio'}->next_seq();
  return @list;
}

sub PRINT {
  my $self = shift;
  $self->{'seqio'}->write_seq(@_);
}

1;

