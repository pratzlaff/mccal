package MCCal::Misc;

use 5.006;
use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw(
		     discrete_draw
		     gauss_truncated
		     gauss_truncated_one
		     hipd_interval
		     modalpoint
		     _interpolate
		     parse_opts
		     perturb_cspline
		     perturb_none
		     perturb_parabola
		     read_specfile
		     read_twocol_bin
		  );

use Carp;
use Log::Any '$log';
use PDL;

=head1 NAME

MCCal::Misc

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SUBROUTINES/METHODS

=cut


# FIXME: not returning nodes
sub perturb_parabola {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  # FIXME: haven't updated to enforce $maxdiff
  my ($e, $emin, $eminvar, $emax, $emaxvar, $maxdiff, $n, $opts) = @_;
  my %opts = %$opts;

  if ($opts{plog}) {
    $_ = log($_) for $e, $emin, $emax;
  }

  my $pert = ones($e->type, $e->nelem);
#  my $pert = ones($e->type, $e->nelem, $n);

  # simple parabolas for now

  # vertex energy has uniform distribution from (emin, emax)
  my $vertex = $rng->get_uniform($n) * ($emax-$emin) + $emin;

  my $dev = $eminvar + ($emaxvar-$eminvar)/($emax-$emin)*($vertex-$emin);

  my $offset;
  # vertex deviation has uniform distribution from (-dev, dev)
  for ($opts{vertexoffset}) {
    $_ eq 'uniform' and $offset = $rng->get_uniform($n) * $dev * 2 - $dev, last;
    $_ eq 'normal' and $offset = gauss_truncated($n, $dev, float), last;
    croak "unrecognized vertex offset distribution = '$_'";
  }

  # calculate coeffs required for maximum perturbation to be attained at
  # "corners"

  my ($i, $tmp);

  my $coeff_ll = (-$dev - $offset) / ($emin - $vertex)**2;
  my $coeff_lr = (-$dev - $offset) / ($emax - $vertex)**2;
  my $coeffs_low = $coeff_lr->copy;
  $i = which($coeff_ll > $coeff_lr);
  ($tmp = $coeffs_low->index($i)) .= $coeff_ll->index($i);

  my $coeff_ul = ($dev - $offset) / ($emin - $vertex)**2;
  my $coeff_ur = ($dev - $offset) / ($emax - $vertex)**2;
  my $coeffs_high = $coeff_ur->copy;
  $i = which($coeff_ul < $coeff_ur);
  ($tmp = $coeffs_high->index($i)) .= $coeff_ul->index($i);

  # now choose coefficients with uniform distribution
  my $coeffs = $rng->get_uniform($n) * ($coeffs_high-$coeffs_low) + $coeffs_low;

  $pert +=
    $offset->dummy(0) +
      $coeffs->dummy(0) *
	($e->dummy(-1, $n) - $vertex->dummy(0))**2;

  my ($x, $y) = ($vertex, $offset);
  print($x, $y, "\n");
  return $pert, $x, $y;
}

sub perturb_cspline {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($e, $emin, $eminvar, $emax, $emaxvar, $maxdiff, $n, $opts) = @_;
  my %opts = %{$opts};

  if ($opts{plog}) {
    $_ = log($_) for $e, $emin, $emax;
  }

  my $pert = ones($e);

  my $x = sequence($opts{spoints}) / ($opts{spoints}-1);
  my $dev = $eminvar + $x * ($emaxvar - $eminvar);

  my $y = zeroes($x);

  $y->set(0, gauss_truncated_one($dev->at(0)));
  for my $i (1..$y->nelem-1) {
    do {
      $y->set($i, gauss_truncated_one($dev->at($i)));
      } while (abs($y->at($i)-$y->at($i-1)) > $maxdiff);
  }

  $x *= $emax-$emin;
  $x += $emin;

=begin comment

  my $x = sequence($opts{spoints}) * ($emax-$emin) / ($opts{spoints}-1) + $emin;
  my $dev = $eminvar + ($emaxvar-$eminvar)/($emax-$emin)*($x-$emin);
#  for my $i (0..$n-1) {

  my $y = $rng->get_uniform($opts{spoints}) * 2 * $dev - $dev;

=cut

=begin comment

how feasible is this as a 15 minute implementation:
restrict the spline control points such that they also have to satisfy eg
c2=c1+/-(0.5*sigma1)
c3=c2+/-(0.5*sigma2)

etc, where sigma1 etc are the deviations allowed at control point 1.
this would restrict the swings from one side to the next.
am envisaging just an extra if statement.

=cut

=begin comment

  for my $j (1..$y->nelem-1) {
    my $max = $dev->at($j);
    my $min = -$dev->at($j);
    $max =$y->at($j-1)+1.0*$dev->at($j-1) if $max>$y->at($j-1)+1.0*$dev->at($j-1);
    $min =$y->at($j-1)-1.0*$dev->at($j-1) if $min<$y->at($j-1)-1.0*$dev->at($j-1);
    (my $tmp = $y->slice("$j:$j")) .= $rng->get_uniform(1)*($max-$min)+$min;
  }

=cut

  my $spl = PDL::GSL::INTERP->init('cspline', $x, $y);

  $pert += $spl->eval($e);
#    (my $tmp = $pert->slice(",($i)")) += $spl->eval($e);
#  }
  return $pert, $x, $y;
}

sub perturb_cspline_old {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my ($e, $emin, $eminvar, $emax, $emaxvar, $n, $opts) = @_;
  my %opts = %{$opts};

  if ($opts{plog}) {
    $_ = log($_) for $e, $emin, $emax;
  }

#  my $pert = ones($e->type, $e->nelem, $n);
  my $pert = ones($e->type, $e->nelem);

  my $x = sequence($opts{spoints}) * ($emax-$emin) / ($opts{spoints}-1) + $emin;
  my $dev = $eminvar + ($emaxvar-$eminvar)/($emax-$emin)*($x-$emin);
#  for my $i (0..$n-1) {

  my $y = $rng->get_uniform($opts{spoints}) * 2 * $dev - $dev;

=begin comment

how feasible is this as a 15 minute implementation:
restrict the spline control points such that they also have to satisfy eg
c2=c1+/-(0.5*sigma1)
c3=c2+/-(0.5*sigma2)

etc, where sigma1 etc are the deviations allowed at control point 1.
this would restrict the swings from one side to the next.
am envisaging just an extra if statement.

=cut

  for my $j (1..$y->nelem-1) {
    my $max = $dev->at($j);
    my $min = -$dev->at($j);
    $max =$y->at($j-1)+1.0*$dev->at($j-1) if $max>$y->at($j-1)+1.0*$dev->at($j-1);
    $min =$y->at($j-1)-1.0*$dev->at($j-1) if $min<$y->at($j-1)-1.0*$dev->at($j-1);
    (my $tmp = $y->slice("$j:$j")) .= $rng->get_uniform(1)*($max-$min)+$min;
  }

    my $spl = PDL::GSL::INTERP->init('cspline', $x, $y);

  $pert += $spl->eval($e);
#    (my $tmp = $pert->slice(",($i)")) += $spl->eval($e);
#  }
  return $pert;
}

# FIXME: not returning nodes
sub perturb_none {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($e, $emin, $eminvar, $emax, $emaxvar, $maxdiff, $n, $opts) = @_;

  return ones($e->type, $e->nelem), pdl(0), pdl(1);
}

sub read_specfile {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my $file = shift;
  my %spec;
  open my $spec, '<', $file or croak "could not open '$file': $!";
  local $_;
  my ($i, $key);
  while (<$spec>) {
    chomp($_);
    next if /^\s*#/ or /^\s*$/;
    if (/^[a-z]/i) {
      $key=$_;
      $i=0;
      next;
    }

    my ($emin, $eminvar, $emax, $emaxvar, $maxdiff, $edgediff) = split ' ', $_;
    push @{$spec{$key}{emin}}, $emin;
    push @{$spec{$key}{eminvar}}, $eminvar;
    push @{$spec{$key}{emax}}, $emax;
    push @{$spec{$key}{emaxvar}}, $emaxvar;
    push @{$spec{$key}{maxdiff}}, $maxdiff;
    push @{$spec{$key}{edgediff}}, $edgediff;

    croak $key unless defined $maxdiff;

    if ( $i and
	 $spec{$key}{edgediff}[-2] <= ($spec{$key}{emaxvar}[-2]-$spec{$key}{eminvar}[-1])) {
	 #$spec{$key}{emaxvar}[-2] - $spec{$key}{edgediff}[-2] - $spec{$key}{eminvar}[-1] > -.005+1e-3) {
      carp("possible issue in spec for $key, emaxvar_$i=$spec{$key}{emaxvar}[-2], edgediff_$i=$spec{$key}{edgediff}[-2], eminvar_@{[$i+1]}=$spec{$key}{eminvar}[-1]\n");
    }
    ++$i;
  }
  $spec->close;
  return %spec;
}

sub gauss_truncated_one {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my $sigma = shift;

  my $s = zeroes(1) + $sigma;

  my $r;
  do { $r = $rng->ran_gaussian_var($s); }
    while (which(abs($r)>$sigma)->nelem);
  return $r->at(0);
}

# outputs a truncated gaussian of length $n
sub gauss_truncated {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my ($n, $sigma, $type) = @_;

  my $s = zeroes($type, $n) + $sigma;
  my $rand = pdl($type, []);

  while ($rand->nelem < $n) {
    my $r = $rng->ran_gaussian_var($s);
    $rand = append($rand, $r->index(which(abs($r) <= $sigma)));
  }

  return $rand->slice('0:'.($n-1));
}

=head2 discrete_draw

=for ref

FIXME

=cut

sub discrete_draw {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $rng = $::rng;

  my ($k, $prop, $n) = @_;

  my $ddh = $rng->ran_discrete_preproc($prop);
  return $k->index($rng->ran_discrete($ddh, $n));

  # used to do it this way, but now we just let GSL take care of things
  my $p = $prop / $prop->sum;
  my $cp = $p->cumusumover;

  my $r = $rng->get_uniform($n);

  my $i = sumover($r->dummy(0, $cp->nelem) > $cp->dummy(-1,$r->nelem));

  return $k->index($i);
}

sub _interpolate {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($xi, $xn, $yn) = @_;

  my ($yi, $err) = interpolate($xi, $xn, $yn);

  my ($i, $tmp);

  $i = which($err & ($xi<=$xn->at(0)));
  ($tmp = $yi->index($i)) .= $yn->at(0);

  $i = which($err & ($xi>=$xn->at(-1)));
  ($tmp = $yi->index($i)) .= $yn->at(-1);

  return $yi;
}

=head2 parse_opts

=for ref

Intended to allow easy options parsing for subroutines. Routines calling
C<parse_opts()> should be designed to allow a hash reference of option
names and values as their last argument. C<parse_opts()> accepts as
input an array reference. If the final element of that array is a
hash reference, it modifies the array by popping the hash reference,
and returns a hash of option names and values.

As an example, if one wrote a subroutine which allowed an "color" option,
then one would call C<parse_opts()> as follows:

	sub draw_image {
		my %opts = parse_opts(\@_, 'color');

		if (defined $opts{color}) {
			# do something with the color passed
		}

		# now play around with @_ to your heart's content

		return;
	}

All arguments given to C<parse_opts()> following the first are interpreted
as "allowed" option names. If an option is given which is not in this list,
then C<Carp::confess()> is called.

C<parse_opts()> works entirely in lowercase. All options passed are
converted to lowercase. Additionally, "allowed" options are converted
to lowercase before comparing to the option names in the input hash
reference. The output hash will have only lowercase key names.

=cut

sub parse_opts {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my $args = shift;

  my %opts = ();

  if (ref $args->[-1] eq 'HASH') {

    # copy hash, making keys lowercase as we go
    my $href_in = pop @$args;
    %opts = map { lc($_) => $href_in->{$_} } keys %$href_in;

    # ensure only allowed options are given
    if (@_) {

      # hash for convenient lookups of valid options
      my %allowed;
      @allowed{map lc, @_} = ();

      my $caller = (caller 1)[3] || 'main';
      exists $allowed{$_} or confess($caller."() - option '$_' unknown")
	for keys %opts;
    }
  }

  return %opts;
}

sub hipd_interval {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my %opts = parse_opts(\@_, qw( nosort ));

  my ($arr, $level) = @_;
  $arr = $arr->qsort unless $opts{nosort};

  my ($min, $max) = $arr->minmax;

  my $mode = modalpoint($arr, { nosort => 1 });
  my $i;

  # duplicate IDL variable names now
  my $nf = $arr->nelem;
  my $ff = $arr;
  my $ii = sequence(long, $arr->nelem);
  my $cff = sequence(double, $arr->nelem)/($arr->nelem-1);
  my $cfmode = interpol($mode, $ff, $cff);
  my $cfmin = $cfmode - $level; $cfmin=0 if $cfmin<0;
  my $cfmax=$cfmode+$level; $cfmax=1 if $cfmax>1;
  my $tmp;
  my $ixmin = interpol($cfmin, $cff, $ii)->long->at;
  $ixmin = 0 if $ixmin<0;
  my $go_on = 1;
  my $k = $ixmin;
  $cfmax = $cff->at($k)+$level;
  my $drng0 = $max-$min;
  my ($hpdm, $hpdp);
  while ($go_on) {
    my $ixmax = $k + long($level*$nf+0.5)->at;
    $ixmax = $nf-1 if $ixmax > $nf-1;
    my $drng = abs($ff->at($k) - $ff->at($ixmax));
    if ($drng<$drng0) {
      $hpdm = $ff->at($k);
      $hpdp = $ff->at($ixmax);
      $drng0 = $drng;
    }
    ++$k;
    $cfmax = $cff->at($k) + $level;
    last if $cff->at($k) > $cfmode;
    last if $cfmax >= 1;
  }

  return $hpdm, $hpdp;
}

=head2 modalpoint

=for ref

Returns the mode of distribution defined by set of unbinned array
values.

This algorithm sorts the input array, divides it into two sets around
the midpoint of its range, chooses the set that contains more
elements, and repeats until the range is sufficiently small. The
midpoint of this final range is the mode.

=for usage

	$modalpoint = modalpoint($pdl, { eps => $epsilon, nosort => 1 });

=cut

sub modalpoint {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my %opts = parse_opts(\@_, qw( eps nosort ));
  my $eps = $opts{eps} ? abs($opts{eps}) : 1e-6;
  my $arr = shift;
  die "insufficient elements in input array" if $arr->nelem < 2;
  return scalar $arr->stats if $arr->nelem < 3;
  $arr = $arr->qsort unless $opts{nosort};
  return modalpoint_recurse($arr, $eps);
}

sub modalpoint_recurse {
  my ($arr, $eps) = @_;

  my ($min, $max) = $arr->minmax;

  return 0.5*($min+$max) if $eps > $max-$min or $arr->nelem == 1;

  my $i = which($arr > 0.5*($max+$min));
  return modalpoint_recurse(
		    ( $i->min > $arr->nelem/2 ?
		      $arr->mslice([0, $i->at(0)-1]) :
		      $arr->index($i)
		    ),
		    $eps
		   );
}

sub read_twocol_bin {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my $file = shift;

  open my $obfh, '<', $file or croak "could not open $file: $!";
  my $data = zeroes(float, 2, (-s $file)/2/4);
  my $num = $obfh->read(${$data->get_dataref}, $data->nelem * 4);
  croak $! unless defined $num;
  $num == $data->nelem * 4 or croak "didn't get correct amount of data";
  $data->upd_data;
  $data->bswap4 unless isbigendian();
  return $data->slice('(0),'), $data->slice('(1),');
}

=begin comment

/data/fubar/SCAR/pro/stat/modalpoint.pro

function modalpoint,array,eps=eps,verbose=verbose, _extra=e
;+
;function modalpoint
;	returns the mode of distribution defined by set of unbinned array values
;
;syntax
;	arraymode=modalpoint(array,eps=eps,verbose=verbose)
;
;parameters
;	array	[INPUT; required] array of values for which mode must be found
;
;keywords
;	eps	[INPUT; default=1e-6] a small number
;	verbose	[INPUT] controls chatter
;	_extra	[JUNK] here only to prevent crashing the program
;
;description
;	sort the array, divide into two sets around midpoint of range,
;	choose the set that contains more elements, and repeat until
;	the range is sufficiently small, and declare the midpoint of
;	the range to be the mode
;
;example
;	for i=0,20 do print,modalpoint(randomn(seed,10000L)+i)
;
;history
;	translated to IDL by Vinay Kashyap from C code written for BEHR
;	  by Taeyoung Park c.2003 (MarMMVI)
;-

;	usage
ok='ok' & np=n_params() & na=n_elements(array)
if np eq 0 then ok='Insufficient parameters' else $
 if na eq 0 then ok='Input array is undefined' else $
  if na lt 2 then ok='Array must have at least 2 elements'
if ok ne 'ok' then begin
  print,'Usage: arraymode=modalpoint(array,eps=eps,verbose=verbose)'
  print,'  return mode of array'
  if np ne 0 then message,ok,/informational
  return,!values.F_NAN
endif

;	inputs and some special cases
if na lt 3 then return,mean(array)
ok=where(finite(array) ne 0,mok)
if mok eq 0 then return,!values.F_NAN
arr=array[ok] & os=sort(arr) & arr=arr[os]
;
vv=0L & if keyword_set(verbose) then vv=long(verbose[0])>1L
;
if keyword_set(eps) then epsilon=double(eps[0]) else epsilon=1d-6

;	step through the array and find the mode
go_on=1
narr=n_elements(arr) & amax=max(arr,min=amin,/nan)
while go_on do begin
  if vv gt 10 then print,strtrim(narr,2)+'.. ',format='($,a)'
  o1=where(arr gt 0.5*(amin+amax),mo1)
  if mo1 eq 0 or mo1 eq narr then message,'BUG?'
  if o1[0] gt narr/2 then tmparr=arr[0:o1[0]-1L] else tmparr=arr[o1]
  if vv gt 100 then print,narr/2,mo1,o1[0]
  arr=tmparr
  narr=n_elements(arr) & amax=max(arr,min=amin,/nan)
  if narr eq 1 then go_on=0	;stop when there is only one element
  if amax-amin lt epsilon then go_on=0	;stop when range gets too small
endwhile

return,0.5*(amin+amax)
end

=end comment

=head1 AUTHOR

Peter Ratzlaff, C<< <pratzlaff at cfa.harvarf.edu> >>

=cut

1; # End of MCCal::Misc
