#! /bin/bash

perl=/proj/axaf/bin/perl
[[ $(hostname) =~ legs|milagro ]] && perl=perl

n=5
infile=simple.arf
infile_base=`basename $infile | sed 's/\..*//' `
outdir=./out
specfile=../data/simple.spec
title='Simple Simulated ARFs'
opts='--speconly --perturb=cspline'
popts='--ratio'

rm -rf "$outdir"
mkdir -p "$outdir"

for i in $(seq 1 $n)
do
    outfile=`printf "%s/%s_%03d.arf" $outdir $infile_base $i`
    "$perl" ../bin/arfmod "$specfile" "$infile" "$outfile" $opts
done

"$perl" plot_arfs.pl "$infile" "$outdir" --title "$title" --dev "$imagedir/${infile_base}"_simulated_arfs.png/png $popts
