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
the `--speconly` and `--specrows` options. An example `specfile` for
XMM is included in the `data` directory, and a run for EPIC pn might be
```
perl bin/arfmod data/xmm.spec epic_pn.arf epic_pn_mod.arf \
     --speconly --specrows=mm,contam,opfm,epicpn
```

The results of 30 runs for this case:
<image src="images/epic_pn_simulated_arfs.png" />
