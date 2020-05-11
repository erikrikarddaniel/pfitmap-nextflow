#!/usr/bin/env sh

rm -r work/ .nextflow.log* ; nextflow run main.nf -with-singularity pfitmap-nextflow.simg -profile test
