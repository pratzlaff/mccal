#! /usr/bin/perl -w
use strict;

=head1 NAME

plot_arfs.pl - Plot ARF simulations

=head1 SYNOPSIS

perl plot_arfs.pl arffile outdir [options]

=head1 DESCRIPTION

Plot ARF simulations created by C<example_arfmod.sh>

=head1 OPTIONS

=over 4

=item --dev=s

PGPLOT output device.

=item --title=s

Plot title.

=item --wav

Plot wavelength instead of energy.

=item --[xlow,xhigh,ylow,yhigh]=f

Specify one or more axis limit.

=item --ylog

Make Y axis logarithmic.

=item --n=i

Plot only the first C<n> simulated ARFS.

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> May 2025

=head1 SEE ALSO

perl(1).

=cut

my $version = '0.1';

use Config;
use Carp;
use File::Basename;
use FindBin;
use PDL;
use PGPLOT;

use Chandra::Tools::Common qw( read_bintbl_cols );

use Getopt::Long;
my %default_opts = (
		    dev => '/xs',
		    n => 30,
		    lw => 3,
		    dlw => 2, # line width for default curve
		    dci => 1, # color index for default curve
		    plw => 1, # line width for perturbed curves
		    pci => 15, # color index for perturbed curves
		    allci => 1, # cycle through color indices
		    ch => 1.5,
		    title => 'Simulated ARFs',
		    );
my %opts = %default_opts;
GetOptions(\%opts,
	   'help!', 'version!', 'debug!',
	   'dev=s', 'n=i', 'lw=f', 'ch=f',
	   'dlw=f', 'dci=i', 'plw=f', 'pci=i', 'allci!',
	   'title=s', 'wav!', 'ratio!',
	   'ylog!', 'ylow=f', 'yhigh=f', 'xlow=f', 'xhigh=f',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV == 2 or die "Usage: $0 arf dir [options]\n\tTry --help for more information.\n";
my ($arf, $dir) = @ARGV;

my $barename = fileparse($arf, qr/\.[^.]*/);

my @arfs = glob("$dir/${barename}_[0-9][0-9][0-9].arf");
my $n = @arfs;
$n = $opts{n} if $opts{n} < $n;

my @ci;
if ($opts{allci}) {
  @ci = (2..15);
  @ci = ( @ci ) x (int($n / @ci) + 1);
}
else { @ci = ($opts{pci}) x $n; }

my ($elo, $ehi, $specresp_) = read_bintbl_cols($arf, 'energ_lo', 'energ_hi', 'specresp', { extname => 'specresp' });
my $e = ($elo+$ehi)/2;

my $specresp = zeros($e->nelem, $n+1);

$specresp->slice(',0') .= $specresp_;

for my $i (1..$n) {
  my ($specresp_) = read_bintbl_cols($arfs[$i-1], 'specresp', { extname => 'specresp' });
  $specresp->slice(",$i") .= $specresp_;
}

my ($x, $y, $axis, $xlabel, $ylabel) = ($e->log10,
					$specresp,
					10,
					'Energy (keV)',
					'EA (cm\u2\d)'
				       );

if ($opts{wav}) {
  $x = 12.398/$e;
  $xlabel = '\gl';
  $axis = 0;
}

if ($opts{ratio}) {
  $y = $specresp / $specresp->slice(',0');
  $ylabel = 'Ratio';
}  

if ($opts{ylog}) {
  $axis += 20;
  $y = $y->log10;
}

my ($xlow, $xhigh) = $x->minmax;
my ($ylow, $yhigh) = $y->minmax;

$ylow *= 0.99;
$yhigh *= 1.01;

$xlow = $opts{xlow} if $opts{xlow};
$xhigh = $opts{xhigh} if $opts{xhigh};

$ylow = $opts{ylow} if $opts{ylow};
$yhigh = $opts{yhigh} if $opts{yhigh};

for my $dev (split(',', $opts{dev})) {
  pgopen($dev) > 0 or die;

  pgslw($opts{lw});
  pgsch($opts{ch});

  pgenv($xlow, $xhigh, $ylow, $yhigh, 0, $axis);

  pglab($xlabel, $ylabel, $opts{title});

  pgslw($opts{plw});

  for my $i (1..$n) {
    pgsci($ci[$i-1]);
    pgline($e->nelem, $x->float->get_dataref, $y->slice(",$i")->float->get_dataref);
  }

  pgslw($opts{dlw});
  pgsci($opts{dci});
  pgline($e->nelem, $x->float->get_dataref, $y->slice(',0')->float->get_dataref);

  pgclos();
}

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}

sub fix_specresp {
  my $y = $_[0]->copy;
  return $y;
  my $limit = 0.1 * $y->max;
  my $i = $y<$limit;
  $y->where($i) .= $limit;
  return $y;
}
