# VMCLitio
# This is a Fortran Code for the Optimization of Explicitely-Correlated Compact Wave Functions for the Study of Lithium Isoelectronic Suquence in the Quartet State (1s2s2p)


clone as:

git clone https://github.com/djuliannader/VMCLitio.git

compile as:

gfortran -o name.exe src/MetropolisLiQuartetPV3_f.f -w

gfortran -o name.exe src/MetropolisLiQuartetPV3_cusp_f.f -w

run as:

name.exe < input.dat > output.dat &

For user manual and data see directory dat
