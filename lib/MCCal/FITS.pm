package MCCal::FITS;

use 5.006;
use strict;
use warnings;

use Exporter 'import';
our @EXPORT_OK  = qw(
		      apply_ratio
		      check_status
		      read_bintbl_cols
		   );

use Astro::FITS::CFITSIO;
use Carp;
use File::Copy;
use Log::Any '$log';
use PDL;

use MCCal::Misc qw( parse_opts );

=head1 NAME

MCCal::FITS - A few routines used by arfmod

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SUBROUTINES/METHODS

=head2 check_status

=for ref

	$retval = check_status($status);

Checks the CFITSIO status variable. If it indicates an error, the
corresponding CFITSIO error message is carp()ed,
and a false value is returned. A true value is returned if the given
status does not indicate an error.

=cut

sub check_status {
  $log->debugf("%s: %s", (caller(0))[3], \@_);
  my $s = shift;
  if ($s != 0) {
    my $txt;
    Astro::FITS::CFITSIO::fits_get_errstatus($s,$txt);
    carp "CFITSIO error: $txt";
    return 0;
  }

  return 1;
}

=head2 read_bintbl_cols

=for ref

	@cols = read_bintbl_cols($file/$fptr,$col1,$col2,...,{options});

Reads the specified columns from a binary table pointed to by the input
filename or C<fitsfilePtr> object, and returns an array of piddles, one
for each column requested. If no columns are specified, all are read
from the binary table. Returns false if an error occured.
Options should be a hash reference of key/value pairs.
For example:

	@cols = read_bintbl_cols('myfile.fits', 'x', 'y',
				 { extname=>'events', status=>1 }
				);

Options available:

=over 4

=item rethash

	# read all columns into a hash
	my %data = read_bintbl_cols('foo.fits[events]', { rethash=>1 });

Specifies that output will be a hash whose keys are the column names read (in
lowercase), with the corresponding piddles for values.

=item colkeys

Specifies that the first output element returned will be a reference to a hash
whose keys are column names, and whose values are themselves hash references which
contain keyword/value pairs for the given column. In addition, an element named
C<_idx> is added specifiying that column's position in the table.

	# read X, Y columns, show colkeys hash
	use Data::Dumper;
	my ($colkeys, %data) = read_bintbl_cols('foo.fits[events]',
						{ rethash=>1, colkeys=>1 }
					       );
	$Data::Dumper::Indent=1;
	print Dumper $colkeys;

	$VAR1 = {
	  'x' => {
	    '_idx' => 11,
	    'ttype' => 'x',
	    'tcuni' => 'deg',
	    'tlmax' => '8.1925000E+03',
	    'tcdlt' => '-1.3666666666667E-04',
	    'tunit' => 'pixel',
	    'tform' => '1E',
	    'tlmin' => '5.0000000E-01',
	    'tcrvl' => '3.2972102733253E+02',
	    'tcrpx' => '4.0965000000000E+03',
	    'tctyp' => 'RA---TAN'
	  },
	  'y' => {
	    '_idx' => 12,
	    'ttype' => 'y',
	    'tcuni' => 'deg',
	    'tlmax' => '8.1925000E+03',
	    'tcdlt' => '1.3666666666667E-04',
	    'tunit' => 'pixel',
	    'tform' => '1E',
	    'tlmin' => '5.0000000E-01',
	    'tcrvl' => '-3.0194870333383E+01',
	    'tcrpx' => '4.0965000000000E+03',
	    'tctyp' => 'DEC--TAN'
	  }
	};

=item extname

Name of the binary table to move to. If a fitsfilePtr object was passed instead
of a filename, the HDU pointer is stored and reset just before C<read_bintbl_cols()>
returns.

=item status

Print to STDERR the amount finished, continually updated, if true.

=item rfilter

CFITSIO-style row filtering specification. Only the rows matching
this filter will be in the output variables.

=item ninc

Number of rows to read incrementally. By default, this number is set
according to C<fits_get_rowsize()> for the table being read.

=item dtypes

C<read_bintbl_cols()> will create the best fit PDL type for each column read.
If the caller wishes, the preferred output datatype for a given column
can be specified using this option. The argument should be a reference to
a hash whose keys are the column names and whose values are C<PDL::Type>
objects of the type wanted. For example:

	($a,$b,$c) = read_bintbl_cols($file,'a','b','c',{dtypes=>{a=>float,c=>short}});

This will force the PDL type of C<$a> to float, and C<$c> to short, while
choosing the best match datatype for C<$b>. It is not possible for the user to
specify dtypes for LOGICAL, ASCII and BIT type columns.

=back

=cut

sub read_bintbl_cols {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my @leave_subs_beg = my @leave_subs_end = ();
  my $cleanup = sub { $_->() for @leave_subs_beg, @leave_subs_end; };
  my $add_clean_beg = sub { push @leave_subs_beg, @_ };
  my $add_clean_end = sub { unshift @leave_subs_end, @_ };

  my $old_packing = Astro::FITS::CFITSIO::PerlyUnpacking(-1);
  Astro::FITS::CFITSIO::PerlyUnpacking(0);
  $add_clean_beg->(sub { Astro::FITS::CFITSIO::PerlyUnpacking($old_packing) }) if $old_packing;

  my %opt = parse_opts(\@_, qw(extname status rfilter ninc dtypes rethash colkeys));

  # first arg is fitsfilePtr or filename
  my $input = shift;
  my $fptr_passed = 0;
  my $fptr;
  my $s = 0;
  if (UNIVERSAL::isa($input,'fitsfilePtr')) {
    $fptr = $input;
    $fptr_passed = 1;
  } else {
    $fptr = Astro::FITS::CFITSIO::open_file($input,Astro::FITS::CFITSIO::READONLY(),$s);
    check_status($s) or
      error("could not open FITS file '$input'"), $cleanup->(), return;
    $add_clean_end->(sub { $fptr->close_file($s=0) } );
  }

  # move to target HDU
  if (exists $opt{extname}){
    my $old_hdu_num;
    $fptr->get_hdu_num($old_hdu_num);
    $fptr->movnam_hdu(Astro::FITS::CFITSIO::BINARY_TBL(),$opt{extname},0,$s);
    check_status($s) or
      error("could not move to binary table '$opt{extname}'"), $cleanup->(), return;
    $add_clean_beg->(sub { $fptr->movabs_hdu($old_hdu_num, undef, $s=0) })
      if $fptr_passed;
  }

  my $type;
  $fptr->get_hdu_type($type, $s);
  $type == Astro::FITS::CFITSIO::BINARY_TBL() or
    error("current HDU is not a binary table extension"), $cleanup->(), return;

  # make data structure describing the columns available
  my $h = $fptr->read_header;
  my $ncols; $fptr->read_key(Astro::FITS::CFITSIO::TINT(),'tfields',$ncols,undef,$s);
  my %cols;
  for (1..$ncols) {
    my ($name, $format);
    $fptr->read_key(Astro::FITS::CFITSIO::TSTRING(), 'ttype'.$_, $name, undef, $s);
    $fptr->read_key(Astro::FITS::CFITSIO::TSTRING(), 'tform'.$_, $format, undef, $s);
    $name = lc $name;
    my ($repeat, $btype) = $format =~ /^(\d*)([A-Z])/ or
      error("do not understand format '$format'"), $cleanup->(), return;
    $repeat = 1 unless length $repeat;
    @{$cols{$name}}{qw( n format repeat btype )} = ($_, $format, $repeat, $btype);
  }

  # read all columns (ordered by their column number), unless argument was given
  my @reqcols = map { lc } @_ ? @_ : sort { $cols{$a}{n} <=> $cols{$b}{n} } keys %cols;


  # construct the colkeys output
  my %colkeys;
  if ($opt{colkeys}) {
    my $hdr = $fptr->read_header;
    $hdr->{$_} =~ s/^'//, $hdr->{$_}=~s/\s*'$// for keys %$hdr; # clean string values
    for my $c (@reqcols) {
      %{$colkeys{$c}} = map { (lc($_)=~/(.*?)\d+$/), $hdr->{$_} }
	grep /^t\D+$cols{$c}{n}$/i, keys %$hdr;
      $colkeys{$c}{_idx} = $cols{$c}{n};
    }
  }

  my %typemap = (
		 'A' => { 'pdl' => undef, 'null_val' => '', },
		 'I' => { 'pdl' => short, 'null_val' => 0, },
		 'J' => { 'pdl' => long, 'null_val' => 0, },
		 'E' => { 'pdl' => float, 'null_val' => 0, },
		 'D' => { 'pdl' => double, 'null_val' => 0, },
		 'X' => { 'pdl' => byte, 'null_val' => 0, },
		 'L' => { 'pdl' => byte, 'null_val' => 0, },
		 'B' => { 'pdl' => byte, 'null_val' => 0, },
		);

  # user cannot specify output types for logical, string and bit columns
  # type specified must be of a PDL::Type token
  croak "dtypes argument must be a hash reference"
    if ($opt{dtypes} and ref $opt{dtypes} ne 'HASH');
  my %user_types = map { lc($_) => $opt{dtypes}{$_} } keys %{$opt{dtypes}};
  for (keys %user_types) {
    if ( $cols{$_}{btype} =~ /^LXA$/i or
	 !UNIVERSAL::isa($user_types{$_},'PDL::Type')
       ) {
      carp("user-specified type for column '$_' being ignored");
      delete $user_types{$_};
    }
  }

  # select datatypes
  for (@reqcols) {
    exists $cols{$_} or
      croak "requested column '$_' not in file";

    exists $typemap{$cols{$_}{btype}} or
      croak "cannot read column type '$cols{$_}{btype}'";

    # the type of piddle we'll read into
    $cols{$_}{ptype} = exists $user_types{$_} ? $user_types{$_} :
      $typemap{$cols{$_}{btype}}{'pdl'};

    # the type we'll tell cfitsio we're reading
    for my $btype ($cols{$_}{btype}) {
      $btype eq 'L' and $cols{$_}{ctype} = Astro::FITS::CFITSIO::TLOGICAL(), last;
      $btype eq 'X' and $cols{$_}{ctype} = Astro::FITS::CFITSIO::TBIT(), last;
      $btype eq 'A' and $cols{$_}{ctype} = undef, last;
      $cols{$_}{ctype} = match_datatype($cols{$_}{ptype});
    }
  }

  my $nrows;
  $fptr->get_num_rows($nrows,$s);
  $s = 0;			# what could possibly go wrong?
  if (!$nrows) {
    $cleanup->();
    my @return = $opt{rethash} ?
      map { $_ => $cols{$_}{btype} eq 'A' ? [] : pdl($cols{$_}{ptype}, []) } @reqcols
      :
      map { $cols{$_}{btype} eq 'A' ? [] : pdl($cols{$_}{ptype}, []) } @reqcols;
    unshift @return, \%colkeys if $opt{colkeys};
    return @return;
  }

  my $ninc;
  if (!$opt{ninc}) { $fptr->get_rowsize($ninc, $s) }
  else { $ninc = $opt{ninc} }
  $ninc = $nrows if $nrows < $ninc;

  # create piddles
  for (@reqcols) {
    if ($cols{$_}{btype} ne 'A') {
      for my $r ($cols{$_}{repeat}) {
	$r == 1 and
	  @{$cols{$_}}{qw( pdl tmppdl )} = (
					    zeroes($cols{$_}{ptype}, $nrows),
					    zeroes($cols{$_}{ptype}, $ninc)
					    ), last;
	$r > 1 and
	  @{$cols{$_}}{qw( pdl tmppdl )} = (
					    zeroes($cols{$_}{ptype}, $r, $nrows),
					    zeroes($cols{$_}{ptype}, $r, $ninc)
					    ), last;
	@{$cols{$_}}{qw( pdl tmppdl )} = undef; # repeat of zero is allowed by spec
      }
    }
    else {
      $cols{$_}{'pdl'} = [];
    }
  }

  # create masks if we'll be row filtering
  my ($good_mask, $tmp_good_mask, $ngood);
  if ($opt{rfilter}) {
    $good_mask = ones(byte,$nrows);
    $tmp_good_mask = ones(byte,$ninc);
    $ngood = 0;
  }

  if ($opt{status}) {
    print STDERR '    ';
    flush STDERR;
  }
  my $rows_done = 0;
  my $pct_done = -1;
  while ($rows_done < $nrows) {


    if ($opt{status}) {
      my $pct_done_update = int(100*$rows_done/$nrows);
      if ($pct_done_update != $pct_done) {
        $pct_done = $pct_done_update;
        printf STDERR "\b\b\b\b %2d%%" , $pct_done;
	flush STDERR;
      }
    }

    my $rows_this_time = $nrows - $rows_done;
    $rows_this_time = $ninc if $rows_this_time > $ninc;

    # row filter
    if ($opt{rfilter}) {
      my $tmp_ngood = 0;
      $fptr->find_rows($opt{rfilter},$rows_done+1,$rows_this_time,$tmp_ngood,${$tmp_good_mask->get_dataref},$s);
      $tmp_good_mask->upd_data;

      check_status($s) or
	error("error filtering rows: rfilter = '$opt{rfilter}'"),
	  $cleanup->(), return;

      (my $t = $good_mask->slice($rows_done.':'.($rows_done+$rows_this_time-1))) .=
	$tmp_good_mask->slice('0:'.($rows_this_time-1));

      $ngood += $tmp_ngood;

      $tmp_ngood > 0 or
	$rows_done += $rows_this_time,
	  next;
    }

    for (@reqcols) {
      next unless $cols{$_}{repeat};

      if ($cols{$_}{btype} ne 'A') {
	$fptr->read_col( $cols{$_}{ctype},
			 $cols{$_}{n},
			 $rows_done+1, 1,
			 $cols{$_}{repeat} * $rows_this_time,
			 0,
			 ${$cols{$_}{tmppdl}->get_dataref},
			 undef,
			 $s);
	$cols{$_}{tmppdl}->upd_data;

	my $slice1 = ($cols{$_}{repeat} > 1 ? ':,' : '').$rows_done.':'.($rows_done+$rows_this_time-1);
	my $slice2 = ($cols{$_}{repeat} > 1 ? ':,' : '').'0:'.($rows_this_time-1);
	(my $t = $cols{$_}{'pdl'}->slice($slice1)) .=
	  $cols{$_}{tmppdl}->slice($slice2);

      } else {			# string type
	my $tmp = [];
	$fptr->read_col(Astro::FITS::CFITSIO::TSTRING(),
			$cols{$_}{n},
			$rows_done+1, 1,
			$rows_this_time,
			0,
			$tmp,
			undef,
			$s,
		       );
	push @{$cols{$_}{'pdl'}}, @$tmp;
      }

      check_status($s) or
	error("error reading FITS data"), $cleanup->(), return;

    }
    $rows_done += $rows_this_time;
  }

  if ($opt{rfilter}) {
    my $good_index = which($good_mask);

    for (@reqcols) {
      next unless $cols{$_}{repeat};

      if ($cols{$_}{btype} ne 'A') {
	my $r = $cols{$_}{repeat};
	if ($r > 1) {
	  my $index = ($good_index->dummy(0,$r)*long($r) + sequence(long,$r))->clump(-1);
	  $cols{$_}{'pdl'} = $cols{$_}{'pdl'}
	    ->clump(-1)
	      ->index($index)
		->reshape($r,$good_index->nelem);
	} else {
	  $cols{$_}{'pdl'} = $cols{$_}{'pdl'}->index($good_index);
	}
      } else {			# string type
	@{$cols{$_}{'pdl'}} = @{$cols{$_}{'pdl'}}[$good_index->list];
      }
    }
  }

  if ($opt{status}) {
    print STDERR "\b\b\b\b100%";
    flush STDERR;
  }

  $cleanup->();

  my @retval = $opt{rethash} ?
    map { $_ => $cols{$_}{'pdl'} } @reqcols :
    map { $cols{$_}{'pdl'} } @reqcols;

  unshift @retval, \%colkeys if $opt{colkeys};

  return @retval;

}

=head2 match_datatype

=for ref

	$cfitsio_type = match_datatype($piddle);
	$cfitsio_type = match_datatype(long); # or short, or float, etc.

PDL datatypes are always guaranteed to be the same size on all architectures,
whereas CFITSIO datatypes (TLONG, for example), will vary on some
architectures since they correspond to the C datatypes on that system. This
poses a problem for Perl scripts which wish to read FITS data into piddles, and
do so in a manner portable to 64-bit architectures, for example.
This routine takes a PDL object or PDL::Types token (returned by float() and friends
when given no arguments), and returns the same-sized CFITSIO datatype, suitable
for passing to routines such as C<fits_read_col()>.

=cut

sub match_datatype {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my $arg = shift;

  my $pdl_type;
  if (UNIVERSAL::isa($arg,'PDL')) {
    $pdl_type = $arg->get_datatype;
  } elsif (UNIVERSAL::isa($arg,'PDL::Type')) {
    $pdl_type = $arg->[0];
  } else {
    croak "argument should be a PDL object or PDL::Type token";
  }

  my $pdl_size = PDL::Core::howbig($pdl_type);

  my @cfitsio_possible_types;
  # test for real datatypes
  if ($pdl_type == float(1)->get_datatype or
      $pdl_type == double(1)->get_datatype
     ) {
    @cfitsio_possible_types = (
			       Astro::FITS::CFITSIO::TDOUBLE(),
			       Astro::FITS::CFITSIO::TFLOAT(),
			      );
  } elsif ($pdl_type == short(1)->get_datatype or
	   $pdl_type == long(1)->get_datatype
	  ) {
    @cfitsio_possible_types = (
			       Astro::FITS::CFITSIO::TSHORT(),
			       Astro::FITS::CFITSIO::TINT(),
			       Astro::FITS::CFITSIO::TLONG(),
			      );
  } elsif ($pdl_type == ushort(1)->get_datatype or
	   $pdl_type == byte(1)->get_datatype
	  ) {
    @cfitsio_possible_types = (
			       Astro::FITS::CFITSIO::TBYTE(),
			       Astro::FITS::CFITSIO::TUSHORT(),
			       Astro::FITS::CFITSIO::TUINT(),
			       Astro::FITS::CFITSIO::TULONG(),
			      );
  } else {
    croak "cannot handle PDL type $pdl_type";
  }


  for (@cfitsio_possible_types) {
    return $_ if $pdl_size == Astro::FITS::CFITSIO::sizeof_datatype($_);
  }

  croak "no CFITSIO type for PDL type $pdl_type";
}

# just print an error message to STDOUT
sub error {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my $txt = shift;
  my $mname = (caller 1)[3];
  defined $mname or $mname = 'main';
  print STDERR "$mname() - $txt\n";
}

sub apply_ratio {
  $log->debugf("%s: %s", (caller(0))[3], \@_);

  my ($infile, $outfile, $extname, $ratio, $history) = @_;

  copy($infile, $outfile) or croak "could not copy '$infile' -> '$outfile': $!\n";

  my $status = 0;
  my $outfptr = Astro::FITS::CFITSIO::open_file($outfile, Astro::FITS::CFITSIO::READWRITE(), $status);
  check_status($status) or croak "error opening output file '$outfile'\n";

  # move to specresp hdu
  $outfptr->movnam_hdu(Astro::FITS::CFITSIO::BINARY_TBL(), $extname, 0, $status);
  check_status($status) or croak "could not move to '$extname' HDU in $outfile\n";

  my %cols = (
	      specresp => { ctype => Astro::FITS::CFITSIO::TDOUBLE(), ptype => double, },
	     );

  for (keys %cols) {
    $cols{$_}{colnum} = undef;
    $outfptr->get_colnum(Astro::FITS::CFITSIO::CASEINSEN(), $_, $cols{$_}{colnum}, $status);
    check_status($status) or croak "'$_' column not found in $extname HDU from $outfile\n";
  }

  # add a HISTORY keyword
  $outfptr->write_history($history, $status);
  check_status($status) or croak "error writing HISTORY entry to $outfile\n";

  my $nrows;
  $outfptr->get_num_rows($nrows, $status);

  for (keys %cols) {
    $cols{$_}{piddle} = zeroes($cols{$_}{ptype}, $nrows);
  }

  $outfptr->perlyunpacking(0);

  for (keys %cols) {
    $outfptr->read_col($cols{$_}{ctype}, $cols{$_}{colnum}, 1, 1, $nrows, 0, ${$cols{$_}{piddle}->get_dataref}, undef, $status);
    $cols{$_}{piddle}->upd_data;
  }
  check_status($status) or croak "error reading data\n";

  # apply ratio, rewrite specresp column
  for (qw( specresp )) {
    (my $tmp = $cols{$_}{piddle}) *= $ratio;
    $outfptr->write_col($cols{$_}{ctype}, $cols{$_}{colnum}, 1, 1, $nrows, $cols{$_}{piddle}->get_dataref, $status);
    check_status($status) or croak "error writing data\n";
  }

  $outfptr->write_chksum($status);
  check_status($status) or croak "error updating checksum in $outfile\n";

  $outfptr->close_file($status);
  check_status($status) or croak "error closing $outfile\n";
}

=head1 AUTHOR

Peter Ratzlaff, C<< <pratzlaff at cfa.harvarf.edu> >>

=cut

1; # End of MCCal::FITS
