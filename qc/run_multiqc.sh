#!/bin/bash

paths=$1
out=$2

multiqc $paths -o $out

