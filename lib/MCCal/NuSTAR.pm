package MCCal::NuSTAR;
use strict;
use warnings;
use feature 'state';

use Exporter 'import';
our @EXPORT_OK = qw( read_spline_bin
		     component_spline
		     read_detabs
		     read_vignet
		     read_arf
		     nuabs
		     perturb_detabs
		     perturb_vignet
		     perturb_effarea
		  );


use FindBin;
use Log::Any '$log';
use PDL;
use PDL::IO::FlexRaw;

use MCCal::FITS 'read_bintbl_cols';

use MCCal;

my $datadir = $MCCal::DATADIR . '/nustar';
my $splines_bin_dir = $datadir . '/splines/bin';

sub read_spline_bin {
  my $subroutine = (caller(0))[3];
  $log->debugf("%s: %s", $subroutine, \@_);

  my $s = shift;

  state %data;

  if (not exists $data{$s}) {
    my $fname = $splines_bin_dir . '/' . $s . 'XSSpline.bin';
    $log->infof("%s: reading splines file '%s'", $subroutine, $fname);
    my ($fX, $fY, $fB, $fC, $fD) = readflex($fname) or die $!;;
    $data{$s} = {};
    @{$data{$s}}{qw/ fX fY fB fC fD /} = ($fX, $fY, $fB, $fC, $fD);
  }

  return @{$data{$s}}{qw/ fX fY fB fC fD /};
}

sub component_spline {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  return component_spline_orig(@_);
  return component_spline_interp(@_);
}

sub component_spline_interp {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my ($s, $x) = @_;

  my ($fX, $fY, $fB, $fC, $fD) = read_spline_bin($s);

  return interpol($x, $fX, $fY);
}

sub component_spline_orig {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($s, $x) = @_;

  my ($fX, $fY, $fB, $fC, $fD) = read_spline_bin($s);
  my $fNp = $fX->nelem;

  my ($fXmin, $fXmax) = $fX->minmax;

  my $y = $x->zeroes;

  for my $i (0..$x->nelem-1) {
    my $x = $x->at($i);
    my $fKstep = 0;
    my $fDelta = -1;

    my $klow = 0;

    if ($x <= $fXmin) { $klow = 0; }
    elsif ($x >= $fXmax) { $klow = $fNp-1 }
    else {
      if ($fKstep) {
	# Equidistant knots, use histogramming
	$klow = int(($x-$fXmin)/$fDelta);
	if ($klow < $fNp-1) { $klow = $fNp-1; }
      }
      else {
	my $khig=$fNp-1;
	my $khalf;
	# Non equidistant knots, binary search
	while($khig-$klow>1) {
	  $khalf = int(($klow+$khig)/2);
	  if ($x > $fX->at($khalf)) { $klow=$khalf; }
	  else { $khig=$khalf; }
	}
      }
    }

    # Evaluate now
    my $dx=$x-$fX->at($klow);
    $y->set($i, $fY->at($klow) + $dx *
	    ($fB->at($klow) + $dx * ($fC->at($klow)+$dx*$fD->at($klow) )));
  }
  return $y;
}

sub perturb_effarea {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my $detabs = read_detabs();
  my $vignet = (read_vignet())[2]->slice(',(6)');

  my ($e, $arf) = read_arf(@_);

  $arf *= $detabs * $vignet;
  my $arfperturb = $arf * perturb_vignet($e) * perturb_detabs($e);

  return $e, $arf, $arfperturb;

}

sub gauss_truncated_one {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my $sigma = shift;

  my $r;
  do { $r = $rng->ran_gaussian_var($sigma)->at(); } while abs($r)>$sigma;
  return $r;
}

sub perturb_vignet {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my $energies = shift;

  my ($e_v, undef, $v) = read_vignet();

  my $oaa = 1; # 1 arcmin offaxis, nominal pointing position

  my $v_default = $v->slice(',(6)')->copy; # 10" resolution

  my $perturb = gauss_truncated_one(0.5); # uncertainty is 30"
  my $oaa_perturb = $oaa + $perturb;
  my $ioaaminus = rint($oaa_perturb*6);
  my $ioaaplus = $ioaaminus + 1;

  my $v_minus = $v->slice(",($ioaaminus)")->copy;
  my $v_plus = $v->slice(",($ioaaplus)")->copy;

  my $v_perturbed = $v_minus + ($v_plus - $v_minus) * ($oaa_perturb*6 - $ioaaminus);

  $_ = interpol($energies, $e_v, $_) for $v_default, $v_perturbed;

  my $ratio = $v_perturbed / $v_default;
  my $nani = which($ratio != $ratio);
  if ($nani->nelem) {
    my $msg = sprintf 'perturb_vignet() - found %d / %d NaNs', $nani->nelem, $ratio->nelem;
#    warn $msg;
    (my $tmp = $ratio->index($nani)) .= 0;
  }
  return $ratio;
}

sub perturb_detabs {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my $energies = shift;

  # pt_range=[0.05,0.15]
  # czt_range=[0.2,0.3]

  my ($pt_thick, $czt_thick) = (0.1, 0.25);
  my $sigma = 0.05;

  my $params = [$pt_thick, $czt_thick, 0, 0.9];

  my $default_nuabs = nuabs($energies, $params);

  $pt_thick += gauss_truncated_one($sigma);
  $czt_thick += gauss_truncated_one($sigma);

  $params = [$pt_thick, $czt_thick, 0, 0.9];

  my $ratio = nuabs($energies, $params) / $default_nuabs;
  my $nani = which($ratio != $ratio);
  if ($nani->nelem) {
    my $msg = sprintf 'perturb_detabs() - found %d / %d NaNs', $nani->nelem, $ratio->nelem;
    warn $msg;
    (my $tmp = $ratio->index($nani)) .= 0;
  }
  return $ratio;
}

sub nuabs_old {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($energies, $params) = @_;
  $params = [0.1, 0.25, 0, 0.9] unless defined $params;

  my $flux = $energies->zeroes;

  my ($PtThickness, $CZTThickness, $ZnThickness, $CdRatio) = @$params; # um
  $_ *= 1e-4 for $PtThickness, $CZTThickness, $ZnThickness; # cm

  use constant CdMass => 112.41; # g/mol
  use constant ZnMass =>  65.39; # g/mol
  use constant TeMass => 127.60; # g/mol
  my $CZTMass = (CdMass * 0.5 * $CdRatio) + (ZnMass * 0.5 * (1.0 - $CdRatio)) + (TeMass * 0.5); # g/mol
  my $CdMassRatio = (CdMass * 0.5 * $CdRatio)         / $CZTMass;
  my $ZnMassRatio = (ZnMass * 0.5 * (1.0 - $CdRatio)) / $CZTMass;
  my $TeMassRatio = (TeMass * 0.5)                   / $CZTMass;

  for my $i (0..$energies->nelem-1) {
    my $PtTrans  = exp(-component_spline('Pt', $energies->at($i)) * $PtThickness);
    my $ZnTrans  = exp(-component_spline('Zn', $energies->at($i)) * $ZnThickness);

    my $CdOfCZTTrans  = exp(-component_spline('CdOfCZT', $energies->at($i)) * $CZTThickness * $CdMassRatio);
    my $ZnOfCZTTrans  = exp(-component_spline('ZnOfCZT', $energies->at($i)) * $CZTThickness * $ZnMassRatio);
    my $TeOfCZTTrans  = exp(-component_spline('TeOfCZT', $energies->at($i)) * $CZTThickness * $TeMassRatio);

    my $CZTTrans  = $CdOfCZTTrans * $ZnOfCZTTrans * $TeOfCZTTrans;

    $flux->set($i, $PtTrans * $CZTTrans * $ZnTrans);
  }

  return $flux;
}

sub nuabs {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($energies, $params) = @_;
  $params = [0.1, 0.25, 0, 0.9] unless defined $params;

  my $flux = $energies->zeroes;

  my ($PtThickness, $CZTThickness, $ZnThickness, $CdRatio) = @$params; # um
  $_ *= 1e-4 for $PtThickness, $CZTThickness, $ZnThickness; # cm

  use constant CdMass => 112.41; # g/mol
  use constant ZnMass =>  65.39; # g/mol
  use constant TeMass => 127.60; # g/mol
  my $CZTMass = (CdMass * 0.5 * $CdRatio) + (ZnMass * 0.5 * (1.0 - $CdRatio)) + (TeMass * 0.5); # g/mol
  my $CdMassRatio = (CdMass * 0.5 * $CdRatio)         / $CZTMass;
  my $ZnMassRatio = (ZnMass * 0.5 * (1.0 - $CdRatio)) / $CZTMass;
  my $TeMassRatio = (TeMass * 0.5)                   / $CZTMass;


  my $PtTrans  = exp(-component_spline('Pt', $energies) * $PtThickness);
  my $ZnTrans  = exp(-component_spline('Zn', $energies) * $ZnThickness);

  my $CdOfCZTTrans  = exp(-component_spline('CdOfCZT', $energies) * $CZTThickness * $CdMassRatio);
  my $ZnOfCZTTrans  = exp(-component_spline('ZnOfCZT', $energies) * $CZTThickness * $ZnMassRatio);
  my $TeOfCZTTrans  = exp(-component_spline('TeOfCZT', $energies) * $CZTThickness * $TeMassRatio);

  my $CZTTrans  = $CdOfCZTTrans * $ZnOfCZTTrans * $TeOfCZTTrans;

  $flux = $PtTrans * $CZTTrans * $ZnTrans;

  return $flux;
}

sub read_detabs {
  my $subroutine = (caller(0))[3];
  $log->debugf("%s: %s", $subroutine, \@_);
  my $detabs_file = $datadir . '/nuAdetabs20100101v002.fits';
  $detabs_file = shift if @_;
  $log->infof("%s: reading detabs file '%s'", $subroutine, $detabs_file);
  my ($elo, $ehi, $detabs) = read_bintbl_cols($detabs_file.'[1]',
					      qw/ energ_lo energ_hi detabs/) or die;

  return 0.5*($elo+$ehi), $detabs;
}

sub read_arf {
  my $subroutine = (caller(0))[3];
  $log->debugf("%s: %s", $subroutine, \@_);
  my $arf_file = $datadir . '/nuA20100101v006.arf';
  $arf_file = shift if @_;
  $log->infof("%s: reading NuSTAR ARF '%s'", $subroutine, $arf_file);
  my ($elo, $ehi, $specresp) = read_bintbl_cols($arf_file.'[1]',
					      qw/ energ_lo energ_hi specresp/) or die;

  return 0.5*($elo+$ehi), $specresp;
}

sub read_vignet {
  my $subroutine = (caller(0))[3];
  $log->debugf("%s: %s", $subroutine, \@_);
  state %vignets;

  my $vignet_file = $datadir . '/nuAvign20100101v006.fits';
  $vignet_file = shift if @_;

  return @{$vignets{$vignet_file}}{qw/energ theta vignet/}
    if exists $vignets{$vignet_file};

  $log->infof("%s: reading NuSTAR vignette '%s'", $subroutine, $vignet_file);
  my ($elo, $ehi, $t, $v) = read_bintbl_cols($vignet_file.'[1]',
					     qw/energ_lo energ_hi theta vignet/
					     ) or die;

  # just use azimuth angle=0
  $_ = $_->slice(',(0)')->sever for $elo, $ehi, $t, $v;

  my $e = 0.5 * ($elo + $ehi);
  $v->reshape($e->dims, $t->dims);

  $vignets{$vignet_file} = {
			    energ => $e,
			    theta => $t,
			    vignet => $v,
			   };

  return $e, $t, $v;
}

1;
