#!/bin/bash

if [[ -f plotsrc/plot.gplt ]]
then
    gnuplot plotsrc/plot.gplt
fi