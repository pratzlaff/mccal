package MCCal::Chandra;

use 5.006;
use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw(
		     acis_model_random
		     aeff_orbit_corr
		     contam_corr
		     hrma_model_random
		  );

use Carp;
use Log::Any '$log';
use PDL;

use MCCal::FITS qw( read_bintbl_cols );
use MCCal::Misc qw( discrete_draw _interpolate gauss_truncated );

=head1 NAME

MCCal::Chandra

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use MCCal::Chandra;

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

=cut

{
  my ($e, $ea);
  sub hrma_ea {
    $log->debugf("%s: %s", (caller(0))[3], \@_);
    if (!defined $e) {
      my $hrmafile = shift;
      my ($elo, $ehigh);
      ($elo, $ehigh, $ea) = read_bintbl_cols($hrmafile, qw( energ_lo energ_hi effarea ), { extname => 'axaf_axeffa' }) or croak;
      $e = ($ehigh+$elo)/2;
      $_ = $_->flat for $e, $ea;
    }
    return $e, $ea;
  }
}

# choose a HRMA model from the random selection, return both default and
# varied models
#
# runs only once regardless of the number of invocations
{
my ($energy, $ea_default, $ea);
sub hrma_model_random {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  return ($energy, $ea_default, $ea) if defined $ea;

  my ($frac, $opts) = @_;
  my %opts = %$opts;
  $frac < 1 and $frac >=0 or croak $frac;

  my ($hrmafile, $hrmarandomfile) = @opts{qw( hrmafile hrmarandomfile )};

  # FIXME: hardcoded number of energies
  my $ngrid = 1322;

  -f $hrmarandomfile or croak "could not find $hrmarandomfile";
  my $size = -s $hrmarandomfile or croak "zero size for $hrmarandomfile";

  my $niter = $size / (4 * $ngrid) - 1;

  # just double check
  $size == ($niter+1)*4*$ngrid or croak "$size, $niter, $ngrid";

  open my $fh, '<', $hrmarandomfile or croak $!;

  # energies first
  $energy = zeroes(float, $ngrid);
  $fh->sysread(${ $energy->get_dataref }, $ngrid * 4) == $ngrid * 4 or croak;
  $energy->upd_data;
  $energy->bswap4 if isbigendian();

  # seek forward to our chosen EA curve
  my $index = int($frac * $niter);
  $fh->sysseek($index * $ngrid * 4, 1);

  # read the EA curve
  $ea = zeroes(float, $ngrid);
  $fh->sysread(${ $ea->get_dataref }, $ngrid * 4) == $ngrid * 4 or croak;
  $fh->close;
  $ea->upd_data;
  $ea->bswap4 if isbigendian();

  # get the standard EA, regrid
  my $energy_default;
 ($energy_default, $ea_default) = hrma_ea($hrmafile);
  $ea_default = _interpolate($energy, $energy_default, $ea_default);

  return $energy, $ea_default, $ea, $index;
}
}

# choose an ACIS QE model from the random selection, return both default
# and varied models
#
# runs only once regardless of the number of invocations
{
my ($energy, $qe_default, $qe, $index);
sub acis_model_random {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  return ($energy, $qe_default, $qe, $index) if defined $qe;

  my ($frac, $opts) = @_;
  my %opts = %$opts;
  my $acisrandomrdb = $opts{acisrandomrdb};

  $frac < 1 and $frac >=0 or croak $frac;

  my @hdr = MyRDB::rdb_header($acisrandomrdb) or croak;
  my %params;
  for (@hdr) {
    /#\s*(\w+):\s+(.*)/ and $params{$1} = $2;
  }
  for (qw( ngrid niter binfile )) {
    exists $params{$_} or croak "no '$_' header parameter found in $acisrandomrdb";
  }
  my $binfile = $params{binfile};
  my $ngrid = $params{ngrid};
  my $niter = $params{niter};

  # append leading directory path to binfile
  $binfile = $1.$binfile if $acisrandomrdb =~ m!(.*/)[^/]!;

  -f $binfile or croak "could not find $binfile";
  my $size = -s $binfile or croak "zero size for $binfile";

  # consistency check
  $size == ($niter+1)*4*$ngrid or
    croak "size=$size, niter=$niter, ngrid=$ngrid, file integrity test failed";

  open my $fh, '<', $binfile or croak "cannot open $binfile: $!";

  # energies first
  $energy = zeroes(float, $ngrid);
  $fh->sysread(${ $energy->get_dataref }, $ngrid * 4) == $ngrid * 4 or croak;
  $energy->upd_data;
  $energy->bswap4 unless isbigendian();

  # default curve next
  $qe_default = zeroes(float, $ngrid);
  $fh->sysread(${ $qe_default->get_dataref }, $ngrid * 4) == $ngrid * 4 or croak;
  $qe_default->upd_data;
  $qe_default->bswap4 unless isbigendian();

  # seek forward to our chosen QE curve
  $index = int($frac * ($niter-1));
  $fh->sysseek($index * $ngrid * 4, 1);

  # read the EA curve
  $qe = zeroes(float, $ngrid);
  $fh->sysread(${ $qe->get_dataref }, $ngrid * 4) == $ngrid * 4 or croak;
  $qe->upd_data;
  $qe->bswap4 unless isbigendian();

  return $energy, $qe_default, $qe, $index;
}
}

# time is seconds since 1998 (i.e., since mjdref, consistent with
# chandra event time tags)
sub contam_corr {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($time, $pfunc, $opts) = @_;

  # returned uncertainties are fractional
  my ($e, $taus, $uncer) = taus($time, $opts);

  my $sumsub = sub {
    my $sum = shift()->copy;
    $sum += $_ for @_;
    return $sum;
  };

  my $trans_orig = exp(-$sumsub->(@$taus));

  for (0..$#{$taus}) {
    my $corr = gauss_truncated(1, $uncer->[$_], $e->type)->at(0);
    $taus->[$_] *= (1+$corr);
  }

  my $trans_corr = exp(-$sumsub->(@$taus));

  my $corr_factor = $trans_corr / $trans_orig;
  return $e, $corr_factor;

=begin comment

  my ($e, $c, $o, $f, $fluff, $cerr, $oerr, $ferr, $flufferr) = taus($time);

  my $emin = $e->at(0) - 1e-3;
  my $emax = $e->at(-1) + 1e-3;

  # now we create perturbations of each component
  my $cp = $c * $pfunc->($e, $emin, $cerr, $emax, $cerr, 1, $opts);
  my $op = $o * $pfunc->($e, $emin, $oerr, $emax, $oerr, 1, $opts);
  my $fp = $f * $pfunc->($e, $emin, $ferr, $emax, $ferr, 1, $opts);
  my $fluffp = $fluff * $pfunc->($e, $emin, $flufferr, $emax, $flufferr, 1, $opts);

  ($c, $e, $o, $f, $fluff) = tau_tcorr($time, $contamtdep, $c, $e, $o, $f, $fluff);

  my $trans_orig = exp(-($c+$o+$f+$fluff));
  my $trans_corr = exp(-($cp+$op+$fp+$fluffp));

  my $corr_factor = $trans_corr / $trans_orig;
  return $e, $corr_factor;

=cut

}

{
  my ($e, @taus, @frac_err);
  # input time is seconds since 1998 (e.g., chandra event times)
  # returned uncertainties are fractional
  sub taus {
    $log->debugf("%s: %s", (caller(0))[3], \@_);

    my ($time, $opts) = @_;
    my %opts = %$opts;
    my ($contamtdep, $contamtau, $tausigma, $factortausigma) =
      @opts{qw( contamtdep contamtau tausigma factortausigma )};

    # read optical depth table
    if (!defined $e) {
      open my $fh, '<', $contamtau or croak "error opening $contamtau: $!";
      local $_;
      while (<$fh>) {
	# found line containing fractional uncertainties
	last if /^\#\s+\d+\.\d+/;
      }

      # extract uncertainties;
      s/^#\s+//;
      @frac_err = split;

      @frac_err = split(',', $tausigma) if $tausigma;
      @frac_err = ($frac_err[0])x4 if @frac_err==1;

      if ($factortausigma) {
	$_ *= $factortausigma for @frac_err;
      }

      # rest of data are optical depths...
      ($e, @taus) = rcols $fh;
      $fh->close;
    }

    return $e, [ tau_tcorr($time, $contamtdep, @taus) ], [ @frac_err ];
  }
}

{
  my ($time, @tcorr);
  # input time is seconds since 1998 (e.g., chandra event times)
  sub tau_tcorr_factors {
    $log->debugf("%s: %s", (caller(0))[3], \@_);

    my ($t, $contamtdep) = @_;

    # read in the time-dependent factors file if called for first time
    if (! defined $time) {
      ($time, @tcorr) = rcols $contamtdep;
    }

    # return the factors at the requested time
    return map { interpol $t, $time, $_ } @tcorr;
  }

}

sub tau_tcorr {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($time, $contamtdep, @tau) = @_;
  my @tcorr = tau_tcorr_factors($time, $contamtdep);
  @tau == @tcorr or croak scalar(@tau) . ' versus ' . scalar(@tcorr);
  return map { $tau[$_] * $tcorr[$_] } 0..$#tau;
}

=head2 aeff_orbit_cor

=for ref

Choose one of the various HRMA models at random, return ratio
relative to the "f" model.

=cut


sub aeff_orbit_corr {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($e, $opts) = @_;
  my %opts = %{$opts};
  my $aefforbitdir = $opts{aefforbitdir};

  my %prob = (
	      'aeff-orbit-200809-01a.fits' =>  1,
	      'aeff-orbit-200809-01b.fits' =>  1,
	      'aeff-orbit-200809-01c.fits' =>  1,
	      'aeff-orbit-200809-01d.fits' =>  1,
	      'aeff-orbit-200809-01e.fits' =>  1,
	      'aeff-orbit-200809-01f.fits' =>  2,
	      'aeff-orbit-200809-01g.fits' =>  1,
	      'aeff-orbit-200810-01v.fits' =>  0.5,
	      );

  my @f = ('aeff-orbit-200809-01a.fits',
	   'aeff-orbit-200809-01b.fits',
	   'aeff-orbit-200809-01c.fits',
	   'aeff-orbit-200809-01d.fits',
	   'aeff-orbit-200809-01e.fits',
	   'aeff-orbit-200809-01f.fits',
	   'aeff-orbit-200809-01g.fits',
	   'aeff-orbit-200810-01v.fits',
	   );

  my @p = ( 1,
	    1,
	    1,
	    1,
	    1,
	    2,
	    1,
	    0.5,
	    );
  my $findex = discrete_draw(sequence(long,scalar @f), pdl(\@p), 1)->at(0);
  my $file = $f[$findex];
  my $ffile = 'aeff-orbit-200809-01f.fits';

#  carp "aeff_orbit_corr() chose $file";

  return(ones($e), $findex) if $file eq $ffile;

  $_ = $aefforbitdir . '/' . $_ for $file, $ffile;
  for ($file, $ffile) {
    -r $_ or croak "maybe --aefforbitdir=s should be used, cannot read $_";
  }

  my ($elo, $ehi, $ea) = read_bintbl_cols($file, qw( energ_lo energ_hi effarea ), { extname => 'axaf_axeffa'}) or croak "error reading $file";
  $ea = _interpolate($e, (0.5*($elo+$ehi))->slice(',(0)'), $ea->slice(',(0)'));

  my ($elof, $ehif, $eaf) = read_bintbl_cols($ffile, qw( energ_lo energ_hi effarea ), { extname => 'axaf_axeffa'}) or croak "error reading $ffile";
  $eaf = _interpolate($e, (0.5*($elof+$ehif))->slice(',(0)'), $eaf->slice(',(0)'));

  return $ea / $eaf, $findex;
}

=head1 AUTHOR

Peter Ratzlaff, C<< <pratzlaff at cfa.harvarf.edu> >>


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc MCCal::Chandra


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

This software is Copyright (c) 2025 by Peter Ratzlaff.

This is free software, licensed under:

  The Artistic License 2.0 (GPL Compatible)


=cut

1; # End of MCCal::Chandra
