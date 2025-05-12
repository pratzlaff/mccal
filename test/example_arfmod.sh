#! /bin/bash

#
# Using one of the ARFs in ./arfs, this script creates a number of ARF
# simulations in ./out and then overplots them all.
#
# Examples runs are...
# 
# ./example_arfmod.sh simple 5
#
# ./example_arfmod.sh epic_pn 30
#
# ./example_arfmod.sh acis_none 30
#
# ./example_arfmod.sh hrcs_letg 30
#

set -e
set -o pipefail

[ $# -eq 2 ] || {
    echo "Usage: $0 inst n" 1>&2
    exit 1
}

inst="$1"
n="$2"

outdir=./out
rm -rf "$outdir"
mkdir -p "$outdir"

perl=perl
specdir=../data
arfdir=./arfs
imagedir=../images
mkdir -p $imagedir

arffile="$arfdir/${inst}.arf"
[ -f "$arffile" ] || {
    echo "ARF file not found: '$arffile'" 1>&2
    exit 1
}

declare -A specs opts popts titles

specs=(
    [simple]=simple
    [hrcs_letg]=chandra
    [acis_none]=chandra
    [epic_pn]=xmm
)
[ -z "${specs["$inst"]}]" ] && {
    echo "unknown instrument configuration: '$inst'" 1>&2
    exit 1
}

specfile="$specdir/${specs["$inst"]}.spec"
[ -f "$specfile" ] || {
    echo "Specfile not found: '$specfile'" 1>&2
    exit 1
}

opts=(
    [simple]='--speconly'
    [hrcs_letg]=
    [acis_none]=
    [epic_pn]="--speconly -specrows=mm,contam,opfm,epicpn"
)

popts=(
    [simple]='--ratio'
    [hrcs_letg]='--wav'
    [acis_none]=
    [epic_pn]=
)

titles=(
    [simple]='Simulated Simulated ARFs'
    [hrcs_letg]='HRC-S/LETG Simulated ARFs'
    [acis_none]='ACIS/NONE Simulated ARFs'
    [epic_pn]='EPIC pn Simulated ARFs'
)

arffile_base=`basename "$arffile" | sed 's/\..*//' `

for i in $(seq 1 $n)
do
    outfile=`printf "%s/%s_%03d.arf" $outdir $arffile_base $i`
    "$perl" ../bin/arfmod "$specfile" "$arffile" "$outfile" ${opts["$inst"]}
done

pngdev="$imagedir/${arffile_base}"_simulated_arfs.png/png

"$perl" plot_arfs.pl "$arffile" "$outdir" --title "${titles["$inst"]}" ${popts["$inst"]} #--dev "$pngdev"
