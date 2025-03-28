#! /bin/bash

n=30
infile=acis_none.arf
infile_base=`basename $infile | sed 's/\..*//' `
outdir=./out2

rm -rf $outdir
mkdir -p $outdir

for i in $(seq 1 $n)
do
    outfile=`printf "%s/%s_%03d.garf" $outdir $infile_base $i`
    perl ../bin/arfmunge $infile $outfile --speconly
done

