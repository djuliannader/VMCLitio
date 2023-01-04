# VMCLitio
# This is a Fortran Package for the Optimization of Compact Wave Functions for the Litium Isoelectronic Suquence in the Quartet State (1s2s2p)


compile as:

gfortran -o name.exe src/MetropolisLiQuartetPV3_f.f -w

gfortran -o name.exe src/MetropolisLiQuartetPV3_cusp_f.f -w

run as:

name.exe < input.dat > output.dat &