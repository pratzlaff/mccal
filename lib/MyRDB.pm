package MyRDB;

use strict;
use warnings;

use Exporter 'import';

our @EXPORT_OK = qw( max );
our @EXPORT = qw( rdb_cols rdb_col_names rdb_header );
our $VERSION = '1.00';

use Carp;
use Log::Any '$log';

our $GZIP = undef;
our $BZIP2 = undef;

#
# - return names of columns in RDB table
# - works with gz and bz2 files
#
sub rdb_col_names {

  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  my $f=_process_filename(shift);

  if ($f ne '-') {
    open(F,$f) or confess("cannot open file '$f': $!\n");
  } else {
    *F = *STDIN;
  }

  local $_;
  while (<F>) {
    /^\#/ or last;
  }
  close F unless $f eq '-';

  chomp;
  $_ or confess "file '$f' not valid RDB format";
  return split /\t/;
}

sub col_names {
  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  print STDERR <<EOP;
MyRDB::col_names() is deprecated, please use MyRDB::rdb_col_names()
EOP
  return rdb_col_names(@_);
}

#
# Return array of header lines, newlines and all.
#
sub rdb_header($) {
  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  my $f=_process_filename(shift);
  my @out=();

  if ($f ne '-') {
    open(F,$f) or confess "cannot open file '$f': $!";
  } else {
    *F = *STDIN;
  }

  local $_;
  while (<F>) {
    /^\#/ or last;
    push @out, $_;
  }

  close(F) unless $f eq '-';
  return @out;
}

########################################
# Return columns from an RDB table. Input is file and list of column names,
# output is list of references to contents of those columns.
# Does absolutely no consistency checks on the input file.
#
sub rdb_cols {
  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  my $f = _process_filename(shift);
  my @cols=@_;
  local $_;
  my @output=(); map { push @output, [] } @cols; # make output array

  # get indices of the columns 
  my @file_cols=rdb_col_names($f) or
    confess "no column names found in file '$f'";
  my %hpos; @hpos{@file_cols}=(0..$#file_cols);
  my @col_indicies=();
  foreach (@cols) {
    (defined $hpos{$_}) and do { push @col_indicies, $hpos{$_}; next };
    confess "column '$_' not in file '$f'";
  }
  my $split_limit=max(@col_indicies)+2;

  # make the arrays
  if ($f ne '-') {
    open(FILE,$f) or confess "cannot open file '$f': $!";

    # if input is STDIN, then rdb_col_names() already did this for us
    while (<FILE>) {
      /^\#/ or last;
    }				# dispose of header
    <FILE>;
  } else {
    *FILE = *STDIN;
    <FILE>;
  }

  while (<FILE>) {
    chop;
    my @temp=split /\t/,$_,$split_limit;
    foreach my $i (0..$#cols) {
      push @{ $output[$i] }, $temp[$col_indicies[$i]];
    }
  }
  close(FILE) unless $f eq '-';
  return @output;
}

sub max {
  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  my $max=shift;
  for (@_) {
    $max=$_ if $max < $_;
  }
  return $max;
}

sub _process_filename {
  $log->debugf("%s: %s", (caller(0))[3], \@_)
    if $log->is_debug;

  my $f = shift;
  for ($f) {
    /\.gz$/ and do {
      my $gzip = defined $GZIP ? $GZIP : 'gzip';
      $f = "$gzip -dc $f |";
      last;
    };
    /\.bz2$/ and do {
      my $bzip2 = defined $BZIP2 ? $BZIP2 : 'bzip2';
      $f = "$bzip2 -dc $f |";
      last;
    };
  }
  return $f;
}

1;

=head1 NAME

MyRDB - routines for reading RDB files

=head1 SYNOPSIS

use MyRDB;

=head1 DESCRIPTION

Simple routines for reading RDB files.

=head1 ROUTINES

When a routine requires an input file, that file
will automatically be be run through I<gzip> or I<bzip2> if
it has an extension of F<.gz> or F<.bz2>, respectively. In the
event that the location of the above programs is not in
your path, you can set the program location with
C<$MyRDB::GZIP> and C<$MyRDB::BZIP2>.

Using a filename argument of '-' indicates that STDIN should be
used.

=over 4

=item rdb_cols

@col_refs=rdb_cols($file,@cols);

Returns, from RDB file C<$file>, a list of array references, each of
which corresponds to the contents of that element's column name given
in C<@cols>. An example: Read I<TRW_ID> and I<X_Ray_energy> columns
from file F<req_run.rdb>

	# import method names
	use MyRDB qw( rdb_cols );

	# read TRW_IDs and energies
	($trw_id, $energy)=rdb_cols('req_run.rdb', 'TRW_ID', 'X_Ray_energy');

	# print first 10 TRW_IDs and energies
	foreach (0..10) {
    	print "$$trw_id[$_]   $$energy[$_]\n";
	}

=item rdb_col_names

@names=rdb_col_names($file);

Returns list of column names from the given RDB file.

=item rdb_header

@lines=rdb_header($file);

Returns list of RDB header lines (newlines still attached).

=back

=cut
