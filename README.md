mccal
========

The source for this project is written in Perl, and uses the
[PDL](https://pdl.perl.org/) and
[Astro::FITS::CFITSIO](https://metacpan.org/pod/Astro::FITS::CFITSIO)
modules.

This began as an investigation into how Chandra X-ray Observatory
calibration uncertainties affect the resulting effective area, and any
subsequent fitted source model parameters. For the typical Chandra
observation, then, a simple invocation is
```
perl bin/arfmod data/chandra.spec orig.arf mod.arf
```

A number of simulated ACIS/NONE ARFs are plotted below. The solid black
line is the "true" response.
<image src="images/acis_none_simulated_arfs.png" />

And for HRC-S/LETG,
<image src="images/hrcs_letg_simulated_arfs.png" />

Bypassing the Chandra-specific routines, a more general use is possible
with the `--speconly` and `--specrows` options. An example `specfile` for
XMM is included in the `data` directory, and a run for EPIC pn might look like
```
perl bin/arfmod data/xmm.spec epic_pn.arf epic_pn_mod.arf \
     --speconly --specrows=mm,contam,opfm,epicpn
```

The results of 30 runs for this case:
<image src="images/epic_pn_simulated_arfs.png" />

The `specfile` format is a component name (telescope, detector, etc.)
followed by line for estimated calibration uncertainties in each desired energy
band.
```
COMP1
emin1 eminsigma1 emax1 emaxsigma1 maxdiff1 edgemax1
emin2 eminsigma2 emax2 emaxsigma2 maxdiff2 ...

COMP2
...
```
where `emin`, `eminsigma`, `emax`, and `emaxsigma` are
the energy band range and uncertainties at each end of the
band. `maxdiff` is the maximum difference allowed between the
two when distribution samples are drawn. `edgemax` is the maximum
difference allowed between draws for the current `emax` and the next band's
`emin`. Thus, `edgemax` is unnecessary for a component's final listed band.

A very simple example would be

```
OBSERVATORY
0.05 0.20 2.0 0.05 0.10 0.06
2.0 0.10 12.0 0.20 0.07
```

This example is listed in `data/simple.spec`, and ratios of
resulting modified to original responses for five runs are plotted below

<image src="images/simple_simulated_arfs.png" />
