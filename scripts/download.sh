#!/bin/bash

wget --user=1kp-data --password=1kp-rna1 -i wanted_assembly_URLs.txt
bzip2 -dv *.bz2