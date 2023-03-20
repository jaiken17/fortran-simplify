#!/bin/bash

cd build

if command -v ifort &> /dev/null
then
    fortc=ifort
else
    fortc=gfortran
fi

make -s FORT=$fortc