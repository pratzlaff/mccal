package MCCal::Misc;

use 5.006;
use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK = qw(
		     parse_opts
		     hipd_interval
		     modalpoint
		  );

use Log::Any '$log';
use PDL;

=head1 NAME

MCCal::Misc - miscellaneous

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use MCCal::Misc;

    my $foo = MCCal::Misc->new();
    ...

=head1 EXPORT

A list of functions that can be exported.  You can delete this section
if you don't export anything, such as for a purely object-oriented module.

=head1 SUBROUTINES/METHODS

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

=head1 BUGS

Please report any bugs or feature requests to C<bug-mccal-misc at rt.cpan.org>, or through
the web interface at L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=MCCal-Misc>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc MCCal::Misc


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<https://rt.cpan.org/NoAuth/Bugs.html?Dist=MCCal-Misc>

=item * CPAN Ratings

L<https://cpanratings.perl.org/d/MCCal-Misc>

=item * Search CPAN

L<https://metacpan.org/release/MCCal-Misc>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

This software is Copyright (c) 2025 by Peter Ratzlaff.

This is free software, licensed under:

  The Artistic License 2.0 (GPL Compatible)


=cut

1; # End of MCCal::Misc
