# File: "wfi_mccal.inp"
# Version: 2.0
# History: 1.0 - Original version
#          2.0 - Aligned to Athena Calibration Requirement Document v0.8
#                Amended inconsistent assumtion on coating (SiC -> B4C)
#
# Author: M.Guainazzi and the calibration requirement validation team
#
# Content: Input file for the MCCal perturbation method (WFI)
#
# Input: Athena Calibration Requirement document v0.8
#        - CAL-EFF-R-003, Relative effective area on-axis, 3%
#	 - CAL-EFF-R-008, Pixel-to-pixel QE uniformity, 1%
# 	 - CAL-EFF-R-0055, Relative area fine structure, 1%
# 	 - Assumed 10% above 10 keV
#
#	List of edges (for a B4C over-coating on Ir)
# 	- C-K: 0.277 keV
# 	- N-K: 0.409 keV
# 	- O-K: 0.543 keV
#       - Si-K: 1.838 keV
# 	- Ir-M1: 3.1737 keV
# 	- Ir-M2: 2.9087 keV
# 	- Ir-M3: 2.5507 keV
# 	- Ir-M4: 2.1161 keV
# 	- Ir-M5: 2.0404 keV
#
# Notes: contamination not included (4% @0.3 keV)
#
# Legenda
# Col.1: Emin - Minimum energy defining the interval
# Col.2: Emindev - Systematic uncertainty of the effective area at Emin
# Col.3: Emax - Maximum energy defining the interval
# Col.4: Emaxdev - Systematic uncertainty of the effective area at Emax
# Col.5: maxdiff - Maximum difference between minimum and maximum
#  	 	   perturbation (assumed equal to the minimum of Emidnev and
#		   Emaxdev)
# Col.6: Edgeveto - Additional uncertainty at the edge energy
#
# JJD: switched columns 5, 6 to be consistent with code input.
#      Made edgeveto 0.01 for 10 keV rows as it has to be > 0.
#
ALL
0.2 0.04 0.277 0.04 0.04 0.01
0.277 0.04 0.409 0.04 0.04 0.01
0.409 0.04 0.543 0.04 0.05 0.01
0.543 0.04 1.838 0.04 0.04 0.01
1.838 0.04 2.0404 0.04 0.04 0.01
2.0404 0.04 2.1161 0.04 0.04 0.01
2.1161 0.04 2.5507 0.04 0.04 0.01
2.5507 0.04 2.9087 0.04 0.04 0.01
2.9087 0.04 3.1737 0.04 0.04 0.01
3.1737 0.04 10.00 0.04 0.04 0.01
10.0 0.10 12.0 0.10 0.10 0.01
