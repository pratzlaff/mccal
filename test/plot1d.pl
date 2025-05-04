#! /usr/bin/perl -w
use strict;

=head1 NAME

template - A template for Perl programs.

=head1 SYNOPSIS

cp template newprog

=head1 DESCRIPTION

blah blah blah

=head1 OPTIONS

=over 4

=item --help

Show help and exit.

=item --version

Show version and exit.

=back

=head1 AUTHOR

Pete Ratzlaff E<lt>pratzlaff@cfa.harvard.eduE<gt> April 2014

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
		    dlw => 5, # line width for default curve
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
	   'title=s',
	   ) or die "Try --help for more information.\n";
if ($opts{debug}) {
  $SIG{__WARN__} = \&Carp::cluck;
  $SIG{__DIE__} = \&Carp::confess;
}
$opts{help} and _help();
$opts{version} and _version();

@ARGV == 2 or die "usage: $0 arf dir\n";
my ($arf, $dir) = @ARGV;

my $barename = fileparse($arf, qr/\.[^.]*/);

my @arfs = glob("$dir/${barename}_[0-9][0-9][0-9].arf");
my $n = @arfs;
$n = $opts{n} if $opts{n} < $n;

pgopen($opts{dev}) > 0 or die;

pgslw($opts{lw});
pgsch($opts{ch});

my @ci;
if ($opts{allci}) {
  @ci = (2..15);
  @ci = ( @ci ) x (int($n / @ci) + 1);
}
else { @ci = ($opts{pci}) x $n; }

my ($elo, $ehi, $specresp) = read_bintbl_cols($arf, 'energ_lo', 'energ_hi', 'specresp', { extname => 'specresp' });
my $e = ($elo+$ehi)/2;

my ($xlow, $xhigh) = $e->minmax;
my ($ylow, $yhigh) = (0, $specresp->max);
$yhigh*=1.2;

$_ = log10($_) for $xlow, $xhigh;#, $ylow, $yhigh;
pgenv($xlow, $xhigh, $ylow, $yhigh, 0, 10);

pglab('energy (keV)', 'effective area (cm\u2\d)', $opts{title});

pgslw($opts{plw});

for my $i (0..$n-1) {
  my $arf = $arfs[$i];
  pgsci($ci[$i]);
  my ($specresp) = read_bintbl_cols($arf, 'specresp', { extname => 'specresp' });
  pgline($e->nelem, $e->log10->float->get_dataref, $specresp->float->get_dataref);
}

pgslw($opts{dlw});
pgsci($opts{dci});
pgline($e->nelem, $e->log10->float->get_dataref, $specresp->float->get_dataref);

pgclos();

exit 0;

sub _help {
  exec("$Config{installbin}/perldoc", '-F', $FindBin::Bin . '/' . $FindBin::RealScript);
}

sub _version {
  print $version,"\n";
  exit 0;
}
