package XSPEC::Fit;
use strict;

require Exporter;
use vars qw(@EXPORT_OK @EXPORT $VERSION @ISA);
@ISA = qw( Exporter );
@EXPORT_OK = ();
@EXPORT = ();
$VERSION = '1.00';

use File::Path;
use File::Copy;
use File::Temp qw( tempdir );
use Carp;

my $LOGFH;

sub new {
  my $class = (ref $_[0]) ? ref(shift) : shift;

  # the object and its data
  my $me = {

	   #
	   # files and directories
	   #
	   OUTDIR => undef,	# where does the fit output go
	   DATA => undef,	# data file
	   BACK => undef,	# background file
	   RESP => undef,	# response matrix file
	   ARF => undef,	# arf file
	   TEMPLATE => undef,	# input template file...
	   TEMPLATE_TEXT => undef, # ...and/or the contents of such a file
	   BASE => undef,	# base name for output files
	   BASE_DEFAULT => 'xspecfit',
	   OUT_PS => undef,	# Postscript of fit and model
	   OUT_LOG => undef,	# 'show all' log file
	   OUT_IN_XCM => undef,	# input file with all fitting commands
	   OUT_OUT_XCM => undef, # fit output
	   TMPD => undef,	# temporary directory

	   #
	   # some XSPEC fitting params
	   #
	   QUERY => 'no',	# XSPEC query command
	   NIT => 20,		# number of iterations
	   DELTA => 0.01,	# critical change in fit statistic
	   IGNORE => undef,	# XSPEC's ignore parameter
	   ERRORS => [],	# list of lists

	   #
	   # set to change XSPEC params before fitting
	   #
	   FIT_CMDS => '',
	   FINAL_CMDS => '',

	   #
	   # logging parameters
	   #
	   LOGON => 0,
	   LOG_ECHO => 1,
	   LOGFILE => undef,

	  };
  bless $me, $class;
}

# for use in derivative classes
sub _no_insert_template {
  return 0;
}

sub outdir {
  my $me = shift;

  if (@_) {
    my $dir = shift;
    (-d $dir) or mkpath($dir) or
      croak "could not create output directory $dir";
    $me->{OUTDIR} = $dir;
  }
  return $me->{OUTDIR};
}

# update list of all errors to run
sub errors {
  my $me=shift;
  @_ and push @{$me->{ERRORS}}, { N => shift, NAME => shift };
  return $me->{ERRORS};
}

# set data file
sub data {
  return $_[0]->_set_input_files($_[1],'DATA');
}

# set background file
sub back {
  return $_[0]->_set_input_files($_[1],'BACK');
}

# set arf
sub arf {
  return $_[0]->_set_input_files($_[1],'ARF');
}

# set response matrix file
sub resp {
  return $_[0]->_set_input_files($_[1],'RESP');
}

# set template file
sub template {
  return $_[0]->_set_input_files($_[1],'TEMPLATE');
}

# can use contents of a template file instead of or in conjunction with an actual file
sub template_text {
  my $me = shift;
  $me->{TEMPLATE_TEXT} = join '', @_ if @_;
  return $me->{TEMPLATE_TEXT};
}

#
# Without args, returns name of current log file (undef if no logging)
# With arg, sets current log file to that name, after first opening the file
# and making sure any old files are closed.
#
sub logfile {
  my $me = shift;

  if (@_) {
    my $logfile = shift;

    if (! defined $logfile) {	# was asked to close current log file
      close $LOGFH if defined $me->{LOGFILE};
      $me->{LOGFILE} = undef;
      return;
    }

    $me->logfile(undef);
    open(FH,">$logfile") or
      carp("could not open file '$logfile': $!"), return;
    $LOGFH = *FH{IO};
    $me->{LOGFILE} = $logfile;
  }

  return $me->{LOGFILE};
}

#
# write a message to the presently defined logfile and STDERR
#
sub log {
  my ($me, $msg) = @_;

  $me->log_echo and print STDERR $msg;
  defined($me->logfile) and print $LOGFH $msg;
}

#
# control state of logging being echoed to STDOUT
#
sub log_echo {
  my $me = shift;
  $me->{LOG_ECHO} = shift if @_;
  return $me->{LOG_ECHO};
}

# general routine for setting data, background, template, and response files
sub _set_input_files { 
  my ($me, $files, $type) = @_;

  my @files;

  #
  # called with undef for $file arg by some routines
  #
  if (defined $files) {

    for my $file (split ',', $files) {
      # might have to strip off trailing {N} (e.g., data file.fits{4})
      my $file_orig = $file;
      $file =~ s/\{\d+\}$//;

      (-f $file and -s _ and -r _) or
	croak "problems with input file '$file' (type $type)";
      push @files, $file_orig;
    }
    $me->{$type}=join ',', @files;
  }
  return $me->{$type};
}

sub query {
  my $me = shift;
  if (@_) {
    for (lc $_[0]) {
      $_ eq 'no' and $me->{QUERY} = 'no', last;
      $_ eq 'yes' and $me->{QUERY} = 'yes', last;
      croak "valid arguments are 'yes' and 'no' only";
    }
  }
  return $me->{QUERY};
}

sub nit {
  my $me = shift;
  return @_ ? ($me->{NIT} = int(shift)) : $me->{NIT};
}

sub delta {
  my $me = shift;
  return @_ ? ($me->{DELTA} = shift) : $me->{DELTA};
}

sub ignore {
  my $me = shift;
  return @_ ? ($me->{IGNORE} = shift) : $me->{IGNORE};
}

# string of XSPEC commands just before fit
sub fit_cmds {
  my $me = shift;
  return @_ ? ($me->{FIT_CMDS} = join '', @_) : $me->{FIT_CMDS};
}

# string of XSPEC commands just before exiting
sub final_cmds {
  my $me = shift;
  return @_ ? ($me->{FINAL_CMDS} = join '', @_) : $me->{FINAL_CMDS};
}

sub base {
  my $me = shift;
  if (@_) {
    for ($_[0]) {
      m!/! and
	croak "output file basenames cannot have slashes ('/')";
      $me->{BASE}=$_;
    }
  }
  return $me->{BASE};
}

sub fit {
  my $me = shift;

  # make sure data file is specified
#  $me->data or croak "no data file specified";

  # make sure output directory is defined
  $me->outdir or $me->outdir('.');

  # set the output file basename
  if (! $me->base) {
    if ($me->data) { $me->base($me->data =~ m!([^/]*?)(:?\.[^.]*)?$!); } # FIXME
    else { $me->base( $me->{BASE_DEFAULT} ); }
  }

  $me->_mktmpd or
    croak "could not make temporary directory";

  my $cmds='';
  $cmds .= "autosave off\n";
  $cmds .= 'query '.$me->query."\n";

  # insert contents of template file
  if ($me->template && ! $me->_no_insert_template) {
    open(TEMPLATE,$me->template) or
      croak "could not open template ".$me->template.": $!";
    $cmds .= join('',<TEMPLATE>);
    close TEMPLATE;
    chomp $cmds; $cmds .= "\n"; # append a newline
  }
  elsif ($me->template_text) {
    $cmds .= $me->template_text;
  }

  $cmds .= ($me->data) ? 'data '.$me->data."\n" : '';
  $cmds .= ($me->back) ? 'backgrnd '.$me->back."\n" : '';
  $cmds .= ($me->resp) ? 'response '.$me->resp."\n" : '';
  $cmds .= ($me->arf) ? 'arf '.$me->arf."\n" : '';
  $cmds .= ($me->ignore) ? 'ignore '.$me->ignore."\n" : '';

  $cmds .= ($me->fit_cmds) ? $me->fit_cmds : '';
  $cmds .= "renorm\n";
  $cmds .= 'fit '.$me->nit.' '.$me->delta."\n";
  $cmds .= 'save all '.$me->tmpd."/out.xcm\n";
  $cmds .= 'log ' .$me->tmpd."/log.log\n";
  $cmds .= "show all\n";
  $cmds .= "log none\n";
  $cmds .= 'cpd '.$me->tmpd."/ps.ps/ps\n";
  $cmds .= "plot ldata resid\n";
#  $cmds .= 'cpd '.$me->tmpd."/gif.gif/gif\n";
  $cmds .= "plot ldata resid\n";

  # now we do error analysis
  my $i=0;
  for (@{$me->errors}) {
    $cmds .= 'log '.$me->tmpd."/logerr$i.log\n";
    $cmds .= 'error STOPAT '.$me->nit.' '.$me->delta.' MAXIMUM 50.0 '.$_->{N}."\n";
    $cmds .= "log none\n";
    $i++;
  }

    # finish up commands
  $cmds .= ($me->final_cmds) ? $me->final_cmds : '';
  $cmds .= "exit\n";
  $cmds .= "y\n";

  # write XSPEC sessions commands and run XSPEC
  my $fbase=$me->outdir."/".$me->base;
  open(IN,">${fbase}_in.xcm") or
    croak "could not open command file '${fbase}_in.xcm' for writing: $!";
  print IN $cmds;
  close(IN);
  #system("/data/calib1/xanadu/SunOS_5.5_sparc/bin/xspec ${fbase}_in.xcm");
  my $ret = system('xspec', '-', "${fbase}_in.xcm");
  $ret == 0 or croak "error running xspec: returned value = $ret, command was xspec - ${fbase}_in.xcm";

  # clean up temporary files
  move($me->tmpd.'/out.xcm', "${fbase}_out.xcm") or
    croak 'error moving '.$me->tmpd."/out.xcm -> ${fbase}_out.xcm: $!";
  move($me->tmpd.'/log.log', "${fbase}_showall.log") or
    croak 'error moving '.$me->tmpd."/log.log -> ${fbase}_showall.log: $!";
  move($me->tmpd.'/ps.ps', "${fbase}_xspec-fit.ps") or
    croak 'error moving '.$me->tmpd."/ps.ps -> ${fbase}_xspec-fit.ps: $!";
#  move($me->tmpd.'/gif.gif_2', "${fbase}_xspec-fit.gif") or
#    croak 'error moving '.$me->tmpd."/gif.gif_2 -> ${fbase}_xspec-fit.gif";

  # move error files
  $i=0;
  foreach (@{$me->errors}) {
    open(F1,$me->tmpd."/logerr$i.log") or die;
    my @F=<F1>;
    close F1;

    open(F2,"> ${fbase}_".$_->{NAME}."err.log") or die;
    print F2 $F[-2];
    close F2;
    $i++;
  }

#  $me->_rmtmpd;

  return 1;
}

sub tmpd {
  my $me = shift;
  return @_ ? ($me->{TMPD} = shift) : $me->{TMPD};
}

# make a temporary directory
sub _mktmpd {

  my $me = shift;

#  my $tempdir = tempdir( CLEANUP => 1)
#    or die "could not create temporary directory: $!";
  my $tempdir = File::Temp->newdir()
    or die "could not create temporary directory: $!";
  $me->tmpd($tempdir);
  return 1;
}

# remove temporary directory
sub _rmtmpd {
  my $me = shift;
  if (rmtree($me->tmpd)) {
    $me->tmpd(undef);
    return 1;
  } else {
    carp "could not remove temporary directory ".$me->tmpd;
    return;
  }
}

1;  # return true

=head1 NAME

XSPEC::Fit - Class for XSPEC fitting

=head1 SYNOPSIS

	use XSPEC::Fit;

=head1 DESCRIPTION

	$fit = new XSPEC::Fit;

Creates an XSPEC fitting 'object', which can be used to do simple XSPEC fits
in a batch manner. For the XSPEC fit, a data file should be specified,
as well as optional background and response files.

	$fit->data( $datafile );        # data file
	$fit->back( $background_file ); # background file
	$fit->resp( $response_matrix ); # rmf
	$fit->arf( $arf );              # arf

For the XSPEC C<ignore> parameter, use the C<ignore()> method. An example:

	$fit->ignore( "1-20 50-55 350-**" );

The model(s) to fit should be in an XSPEC template file. This template
is usually an XSPEC fit of a similar dataset saved using the XSPEC command

	XSPEC> save mo template

To set the template for the model, use the C<template()> method:

	$fit->template( $XSPEC_template_file );

The contents of the template file are inserted into the list of
input XSPEC commands. Additionally, the C<template_text()> method
can be used to add the contents of such a template file, as opposed
to just giving the filename to C<template()> and letting the module
get the contents itself. Both C<template()> and C<template_text()>
can be used in the same fit, with the C<template()> file contents
being inserted first.

The C<fit()> method does the fit itself:

	$fit->fit;

XSPEC's C<fit> command takes the general form of

	XSPEC> fit niterations delta

where C<niterations> and C<delta> are the number of fit iterations
and critial change in fit statistic respectively. By default, the
number of iterations is 20 and the change in fit statistic is 0.01.
These values can be changed using the C<nit()> and C<delta()> methods.

	$fit->nit( 10 ); # set number of fit iterations to 10
	$fit->delta( .05 ); # set value of critical change in statistic

Output of the fit is in the form of four files:

=over 4

=item XSPEC input commands

A filename of the form F<BASE_in.xcm> containing all commands for the
XSPEC session.

=item XSPEC output fit

Filename of the form F<BASE_out.xcm> containing the fit using the XSPEC
command

	XSPEC> save mo BASE_out

=item XSPEC output log

Filename of form F<BASE_showall.log> containing output of XSPEC command

	XSPEC> show all

=item Postscript plot

Filename of form F<BASE_xspec-fit.ps> containing plot of data and folded
XSPEC model.

=back

The C<BASE> in the above filenames is, by default, taken as the basename
of the input data file (i.e., filename minus extension and leading
directory name). This can be set to a custom value using the C<base()> method:

	$fit->base( $my_very_own_basename );

The above output files are, by default, put into the current working
directory. You can change that directory with the C<outdir()> method:

	$fit->outdir( $my_very_own_output_directory );

Note that not only do the above routines set parameters, they also
return the current value of that parameter.

The user may wish to have additional XSPEC commands put into the
fit. This can be done using the C<fit_cmds()> method.

	$fit->fit_cmds( "thaw 5\n" ); # set XSPEC commands just before fit

These additional commmands are inserted I<just> before the XSPEC fit is
done. Note that it is up to the user to put newlines into commands.
There is also a way of inserting XSPEC commands just before XSPEC is
exited after the fit is done. The necessary method for this is
C<final_cmds()>.

=head1 METHODS

In general, each of the below methods will return the value of a particular
object data item, as well as setting that value. To retrieve the present
value for said data item without changing it, one would usually call the
method without an argument.

=over 4

=item new

The constructor. Simply creates object with data.

=item outdir

Set output directory. Attempts to create the directory if it does not
already exist.

=item query

Change value of XSPEC's C<query> paramater. Valid values are
C<yes> and C<no>. The default is C<no>.

=item data

Set name of XSPEC data file.

=item back

Set name of XSPEC background file.

=item resp

Set name of XSPEC response matrix.

=item template

Set name of an XSPEC input file.

=item template_text

Set contents of an XSPEC input file.

=item base

Set basename of output files.

=item nit

Change number of iterations for fit (default is 20).

=item delta

Change value of critical change in fit statistic (default is 0.01).

=item fit_cmds

Set string of XSPEC commands executed just before the fit is done.

=item final_cmds

Set string of XSPEC commands executed just before exiting XSPEC.

=item fit

Do the XSPEC fit.

=back

=head1 BUGS/FEATURES

FIXME

=head1 AUTHOR

Pete Ratzlaff <pratzlaff@cfa.harvard.edu>

Copyright 2001, Smithsonian Institution

=cut
