import argparse
import astropy.io.fits
import glob
import matplotlib
import matplotlib.pyplot as plt
import os
from pathlib import Path

default_arf = None
default_e = None

def plot_arf(arf, args, symbol='-'):
    global default_arf

    e, specresp = read_arf(arf)
    x = e
    y = specresp
    if args.ratio:
        y = specresp / default_arf
    if args.wav:
        x = 12.398 / e
    plt.plot(x, y, symbol)

def read_arf(arf):
    with astropy.io.fits.open(arf) as hdulist:
        data = hdulist['specresp'].data
        energ_lo = data['energ_lo']
        energ_hi = data['energ_hi']
        specresp = data['specresp']
        return 0.5*(energ_lo+energ_hi), specresp

def get_arfs(arf, dir, n):
    stem = Path(arf).stem
    globstr = dir+'/'+stem+'_[0-9]*.*'
    arfs = sorted(glob.glob(globstr))

    if n:
        arfs = arfs[:n]

    return arfs

def plot_arfs(args):

    global default_arf, default_e

    if args.ratio:
        default_e, default_arf = read_arf(args.arf)

    if args.lw:
        matplotlib.rcParams['lines.linewidth'] = args.lw
        matplotlib.rcParams['axes.linewidth'] = args.lw
    if args.fs:
        matplotlib.rcParams['font.size'] = args.fs

    if args.xlog:
        plt.xscale('log')
    if args.ylog:
        plt.yscale('log')

    arfs = get_arfs(args.arf, args.dir, args.n)
    for arf in arfs:
        plot_arf(arf, args)
    plot_arf(args.arf, args, '-k')

    plt.title(args.title)

    xlabel = 'Energy (keV)'
    if args.wav:
        xlabel = 'λ (Å)'
    ylabel = r'$\mathrm{EA\;(cm^2)}$'
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    if args.xmin:
        plt.xlim(left=args.xmin)
    if args.xmax:
        plt.xlim(right=args.xmax)
    if args.ymin:
        plt.ylim(bottom=args.ymin)
    if args.ymax:
        plt.ylim(top=args.ymax)

    plt.tight_layout()

    if (args.outfile):
        plt.savefig(args.outfile)
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(
        description='Plot ARF and permutations.'
    )
    parser.add_argument('-o', '--outfile', help='Save plot to named file.')
    parser.add_argument('-t', '--title', default='Simulated ARFs', help='Plot title.')
    parser.add_argument('-n', type=int, help='Maximum number of ARFs to plot.')
    parser.add_argument('-w', '--wav', action='store_true', help='Plot vs wavelength.')
    parser.add_argument('-r', '--ratio', action='store_true', help='Plot ratios vs default ARF.')
    parser.add_argument('--xlog', action='store_true')
    parser.add_argument('--ylog', action='store_true')
    parser.add_argument('--xmin', type=float, help='Lower X limit.')
    parser.add_argument('--xmax', type=float, help='Upper X limit.')
    parser.add_argument('--ymin', type=float, help='Lower Y limit.')
    parser.add_argument('--ymax', type=float, help='Upper Y limit.')
    parser.add_argument('--lw', type=float, help='Line widths.')
    parser.add_argument('--fs', type=float, help='Font sizes.')
    parser.add_argument('arf', help='Input ARF.')
    parser.add_argument('dir', help='Directory containing mutated ARFs.')
    args = parser.parse_args()

    plot_arfs(args)


if __name__ == '__main__':
    main()

