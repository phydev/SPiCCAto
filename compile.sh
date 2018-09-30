#!/bin/bash
rm -fr mod
rm run-cells 
mkdir mod
dir='./src'
ifort -O3 -std03 -fast -r8 -g -traceback -module mod $dir/global_m.F90 $dir/misc_m.F90 $dir/derivatives_m.F90 $dir/sim_init_m.F90 $dir/run_cells_m.F90  $dir/main.F90 -o run-cells

#fort -O0 -g -traceback -r8 -module mod $dir/global_m.F90 $dir/misc_m.F90 $dir/derivatives_m.F90 $dir/sim_init_m.F90 $dir/run_cells_m.F90  $dir/main.F90 -o cells3d
