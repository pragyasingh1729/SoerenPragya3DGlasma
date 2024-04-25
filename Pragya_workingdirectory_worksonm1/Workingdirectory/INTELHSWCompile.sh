#!/bin/bash

# MAKES BINARIES THAT WORK FOR CORI HASWELL

module load intel
module load cray-fftw
module load gsl

scp INTELMakefileHSW Makefile

make all
