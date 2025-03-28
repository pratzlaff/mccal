package XSPEC::Output;
use strict;

use vars qw( @EXPORT_OK @EXPORT $VERSION @ISA );

require Exporter;

@ISA = qw( Exporter );
@EXPORT_OK = qw( read_xspec_log read_log );
@EXPORT = qw( );
$VERSION = '1.00';

use Carp;

# read an xspec log file for model jmkmod? made with commands:
#
# XSPEC> log $file
# XSPEC> show all
# XSPEC> log none
#
sub read_xspec_log {

  carp "XSPEC::Output::read_xspec_log() is deprecated, use XSPEC::Output::read_log() instead";

  my @O = read_log($_[0]);

  my $O={};
  for my $href (@O) {
    for my $key (keys %$href)  {
      # save parameter name, value and error
      my ($par,$val)=($href->{$key}{par},$href->{$key}{val});
      my $err = (exists $href->{$key}{err}) ? $href->{$key}{err} : undef;

      if (exists $O->{$par}) {
	$O->{$par}=[ $O->{$par} ] unless
	  (ref $O->{$par}[0]); # not array ref with val/err pair
	push @{$O->{$par}}, [ $val, $err ];
      } else {
	$O->{$par} = [ $val, $err ];
      }
    }
  }

  return $O;
}

sub read_log_xspec11 {

  my $fh = shift;

  my %config;

  while (<$fh>) {
    last if /^\s*Model\s+Fit\s+Model\s+Component\s+Parameter/;

    $config{statistic} = lc($1), next if
      /fit statistic in use.*?(\S*)$/i;

    $config{technique} = $1, next if
      /minimization technique is (\S*)/i;

    $config{criterion} = $1, next if
      /convergence criterion =\s+(\S*)/i;

    @config{qw( telescope instrument type )} = ($1,$2,$3), next if
      /telescope = (\S*) , instrument = (\S*) , channel type = (\S*)/i;

    $config{exptime} = $1, next if
      /file integration time\s+(\S*)/i;

    $config{observed_rate} = $1, next if
      /file observed count rate\s+(\S*)/i;

    $config{predicted_rate} = $1, next if
      /model predicted rate :\s+(\S*)/i;

  }
  (defined $_) or die "file unrecognized";
  <$fh>;				# skip line

  my @output;
  my %errs;		    # keep track of uncertainties by model par
  while (<$fh>) {

    chomp;
    last if /^\s*-------/;
    s/^\s*(.*?)\s*$/$1/; # remove leading and trailing spacings

    # pull out everything before a Unit field
    s/^(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+//;
    length or die "file unrecognized";

    my ($mpar, $fpar, $mcomp, $comp, $par) = ($1,$2,$3,$4,$5);
    @{$output[$mcomp-1]{$par}}{qw( mpar fpar mcomp comp par )} = ($mpar,
								  $fpar,
								  $mcomp,
								  $comp,
								  $par);

    # get value and error
    s/^\S+\s+(.*)$/$1/ unless /^(?:[+-]?)(?=\d|\.\d)\d*(\.\d*)?(?:[Ee](?:[+-]?\d+))?\s/; # remove Unit field, if necessary (regex taken from perldoc -q float, to test for C float)

    /^(\S+)\s+(.*)$/ or die;
    $output[$mcomp-1]{$par}{val} = $1;

    my $err=$2;	# the error portion

    for ($err) {

      # straight uncertainty, save and move on
      m!\+/-\s*(.*)! and
	$errs{$mpar} = $output[$mcomp-1]{$par}{err} = $1, last;

      # parameter was fixed
      /frozen/i and last;

      # dependent
      /=\s*par\s*(\d+)\s*(.*)/i and do {

	my ($dep, $tmp) = ($1, $2);

	last unless exists $errs{$dep}; # dependent param was frozen

	my $factor;
	if ($tmp =~ /\*\s*(.*)/) {
	  $factor = $1;
	}
	else { # additive, or no modification to dependency
	  $factor = 1;
	}

	$errs{$dep} < 0 and $factor = 1; # preserve -1

	$errs{$mpar} = $output[$mcomp-1]{$par}{err} = $errs{$dep} * $factor;

      }, last;

      die $err;
    }
  }

  defined $_ or die "file unrecognized";

  while (<$fh>) {

   @config{ qw( fit_stat nbins )} = ($1, $2), next if
    /\Q$config{statistic}\E\s+=\s+(\S+) using (\d+) PHA bins/i;

   @config{ qw( test_stat )} = ($1), next if
    /test statistic : chi-squared =\s+(\S+) using (\d+) PHA bins/i;

   @config{ qw( reduced_chi_squared dof ) } = ($1, $2), next if
      /reduced chi-squared =\s+(\S+) for\s+(\d+) degrees of freedom/i;

 }

  $output[$_]{CONFIG} = \%config for (0..$#output);
  return @output;
}

sub read_log {
  local $_;
  my $file = shift;

  my $we_opened;
  if (! ref $file) {
    ($file =~ /\.gz$/) and $file = "gzip -dc $file |";
    ($file =~ /\.bz$/) and $file = "bzip2 -dc $file |";
    open(F,$file) or croak "could not open file '$file': $!";
    $we_opened = 1;
  }
  else {
    *F = *{$file}{IO};
  }

  # find a date string, which tells us which xspec version was used,
  # call appropriate parse routine

  my @return;
  while (<F>) {
    if (/^\s+\d+:\d+:\d+\s+\d+-\w+-\d+/) {
      @return = read_log_xspec11(*F{IO});
      last;
    }
    elsif (/^#\w+\s+\w+\s+\d+\s+\d+:\d+:\d+\s+\d+/) {
      @return = read_log_xspec12(*F{IO});
      last;
    }
  }

  # reached end of file without parsing it
  @return or die "contents of '$file' are unrecognized";

  close F if $we_opened;

  return @return;
}

sub read_log_xspec12 {
  my $fh = shift;

  my %config;

  # skip header
  while (<$fh>) {
    last if /^#Model\s+Model\s+Component\s+Parameter/;

    $config{statistic} = lc($1), next if
      /fit statistic in use.*?(\S+)$/i;

    $config{technique} = $1, next if
      /minimization technique:\s+(\S+)/i;

    $config{criterion} = $1, next if
      /convergence criterion =\s+(\S+)/i;

    @config{qw( telescope instrument type )} = ($1,$2,$3), next if
      /telescope:\s+(\S*)\s+instrument:\s+(\S*)\s+channel type:\s+(\S*)/i;

    # assume first exposure time we come across is for data file
    if (/exposure time:\s+(\S*)/i) {
      $config{exptime} = $1 unless exists $config{exptime};
      next;
    }

    $config{counts} = $1, next if
      /source spectrum counts: (\S+)/i;

    $config{predicted_rate} = $1, next if
      /model predicted rate: (\S*)/i;

  }
  (defined $_) or die "file unrecognized";
  <$fh>; # skip line

  my @output;
  my %errs; # keep track of uncertainties by model par
  my $in_data_group;
  while (<$fh>) {

    chomp;
    last if /^#______/;

    # FIXME
    if (/data group: \d+/i) {
      last if $in_data_group;
      $in_data_group = 1;
      next;
    }

    # remove leading and trailing spacings
    s/^#\s*(.*?)\s*$/$1/;

    # pull out everything before a Unit field
    s/^(\d+)\s+(\d+)\s+(\S+)\s+(\S+)\s+//;
    length or die "file unrecognized";

    my ($mpar, $mcomp, $comp, $par) = ($1,$2,$3,$4);
    @{$output[$mcomp-1]{$par}}{qw( mpar mcomp comp par )} = ($mpar,
							     $mcomp,
							     $comp,
							     $par);

    # bool have no unit or uncertainty field
    if (/^(true|false)/) {
      $output[$mcomp-1]{$par}{val} = $1 eq 'true' ? 1 : 0;
      next;
    }

    # remove Unit field, if necessary (regex taken from perldoc -q float, to test for C float)
    s/^\S+\s+(.*)$/$1/ unless /^(?:[+-]?)(?=\d|\.\d)\d*(\.\d*)?(?:[Ee](?:[+-]?\d+))?\s/;

    # get value and error
    /^(\S+)\s+(.*)$/ or die;
    $output[$mcomp-1]{$par}{val} = $1;

    my $err=$2;	# the error portion

    for ($err) {

      # straight uncertainty, save and move on
      m!\+/-\s*(.*)! and
	$errs{$mcomp}{$par} = $output[$mcomp-1]{$par}{err} = $1, last;
#	$errs{$mpar} = $output[$mcomp-1]{$par}{err} = $1, last;

      # parameter was fixed
      /frozen/i and last;

      # FIXME - XSPEC 12.8 just uses, e.g., "= 20", or is it only for tied parameters within the same model?
      /=\s+(\d+)$/ and last;


      # dependent
      /=\s*\w+\[(\d+)\]:(\w+)\s*(.*)/i and do {
 #     /=\s*par\s*(\d+)\s*(.*)/i and do {

	my ($mcomp_dep, $par_dep, $tmp) = ($1, $2, $3);
#	my ($dep, $tmp) = ($1, $2);

	# dependent param was frozen
	last unless exists $errs{$mcomp_dep}{$par_dep};

	my $factor;
	if ($tmp =~ /\*\s*(.*)/) {
	  $factor = $1;
	}
	else { # additive, or no modification to dependency
	  $factor = 1;
	}

	$errs{$mcomp_dep}{$par_dep} < 0 and $factor = 1; # preserve -1
#	$errs{$dep} < 0 and $factor = 1; # preserve -1

	$errs{$mcomp}{$par} =
	  $output[$mcomp-1]{$par}{err} = $errs{$mcomp_dep}{$par_dep} * $factor;
#	$errs{$mpar} = $output[$mcomp-1]{$par}{err} = $errs{$dep} * $factor;

      }, last;

      die $err;
    }
  }

  defined $_ or die "file unrecognized";

  while (<$fh>) {

   @config{ qw( fit_stat nbins )} = ($1, $2), next if
    /fit statistic : \Q$config{statistic}\E\s+=\s+(\S+) using (\d+) PHA bins/i;

   @config{ qw( test_stat )} = ($1), next if
    /test statistic : chi-squared =\s+(\S+) using (\d+) PHA bins/i;

   @config{ qw( reduced_chi_squared dof ) } = ($1, $2), next if
      /reduced chi-squared =\s+(\S+) for\s+(\d+) degrees of freedom/i;
 }

  $output[$_]{CONFIG} = \%config for (0..$#output);
  return @output;
}

1;

=head1 NAME

XSPEC::Output - Module for XSPEC output

=head1 SYNOPSIS

	use XSPEC::Output;

=head1 DESCRIPTION

Tools for playing with XSPEC output. The included routines have been
designed to work with XSPEC v9 and v10 'show all' output.

=head1 ROUTINES

=over 4

=item read_log

	@A=read_log($file);

This is a general purpose routine for reading parameters from an
XSPEC 'show all' log file. The output consists of an array of
hash references, one for each model component. The hash corresponding
to each model component have as keys the parameter names for the model.
The values corresponding to each key are themselves another hash, having
keys:

=over 4

=item mpar

The 'Model par' field in the log file, denoting that parameter's
position in the model.

=item fpar

The 'Fit par' field, vailable only in XSPEC 11 output.

=item mcomp

The 'Model comp' field, denoting that model's position
in the fit.

=item comp

The name of the model (e.g., C<jmkmod4>).

=item par

The name of the parameter in the current model (e.g., C<fano>).

=item val

The value used in the fit.

=back

Additionally, the hash reference will have a C<err> key if the current
parameter was not frozen or coupled with some other parameter. Note
that this value may be -1.0 if the fit failed for some reason, or
the error may even be 0 in other unusual circumstances.

Finally, each hash in the output array will have an additional key,
C<CONFIG>, containing miscellaneous data. An example of these data are:

	$VAR1 = {
	          'nbins' => '444',
	          'criterion' => '9.9999997764826E-03',
	          'observed_rate' => '197.1',
	          'telescope' => 'HXDS',
	          'reduced_fit_stat' => '0.6026659',
	          'fit_stat' => '264.5703',
	          'predicted_rate' => '184.4',
	          'statistic' => 'Chi-Squared',
	          'technique' => 'Lev-Marq',
	          'exptime' => '19.88',
	          'type' => 'PHA',
	          'instrument' => 'fpc_hn'
	        };


=item read_xspec_log

B<WARNING: this is deprecated. You should use C<read_log()>>

	$href=read_xspec_log($file);

Returns a hash reference whose keys are XSPEC model parameter names,
and whose values are array references containing values and
uncertainties (C<undef> if frozen). If a parameter name is not unique
(an example would be addition of two identical models, which
would make each parameter name appear twice), then the
hash value for that parameter name is a reference to an array,
each element of which is its own array of value/error pairs.

To summarize, two scenarios:

	$href=read_xspec_log($file);
	foreach my $key (keys %{$href}) {
		my ($value, $error);
		if (! ref $$href{$key}[0]) { # this is a unique parameter name
			($value, $error) = @{$$href{$key}};
		}
		else {  # non-unique parameter name, use array refs for each instance
			foreach my $i (0..$#{$$href{$key}}) {
				($value, $error) = @{$$href{$key}[$i]};
			}
		}
	}

If this sounds confusing to you (and I know it does) and you have suggestions
for a better interface, then please let me know.

=back

=head1 BUGS/FEATURES

None known, which implies nothing. Please report any strange behaviour (on
the part of this module, that is). Also, I am very interested in hearing
suggestions for additions, improvements, etc.

=head1 AUTHOR

Peter Ratzlaff <pratzlaff@cfa.harvard.edu>

Copyright 2001, Smithsonian Institution

=cut
