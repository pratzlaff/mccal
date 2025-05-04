#! /bin/bash

perl=/proj/axaf/bin/perl
[[ $(hostname) =~ legs|milagro ]] && perl=perl

n=30
infile=acis_none.arf
infile_base=`basename $infile | sed 's/\..*//' `
outdir=./out
specfile=../data/chandra.spec
title='ACIS/NONE Simulated ARFs'
opts="" #--speconly"

\rm -rf "$outdir"
\mkdir -p "$outdir"

for i in $(seq 1 $n)
do
    outfile=`printf "%s/%s_%03d.arf" $outdir $infile_base $i`
    "$perl" ../bin/arfmod "$specfile" "$infile" "$outfile" $opts
done

"$perl" plot1d.pl "$infile" "$outdir" --title "$title" --dev "${infile_base}"_simulated_arfs.png/png
