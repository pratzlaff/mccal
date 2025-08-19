package MCCal;

use 5.006;
use strict;
use warnings;

use Exporter 'import';

our $DATADIR = $INC{'MCCal.pm'};
$DATADIR =~ s!MCCal\.pm$!../data!;

1;
