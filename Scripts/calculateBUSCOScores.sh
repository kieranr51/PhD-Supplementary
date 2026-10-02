#!/bin/bash

# This script goes through and runs eggnog-mapper on all of the proteome files produced by
# the parasitism island pipeline

mkdir -p buscoAssemblyResults
mkdir -p buscoAnnotationResults

for group in parasitism_island_results/*; do
    groupName=$(basename $group)
    for species in $group/*; do
        speciesName=$(basename $species)

        if [ $speciesName != tree.txt ] && [ $speciesName != orthofinder ] && [ ! -d buscoAssemblyResults/$speciesName ]; then
            compleasm run -l eukaryota -t 10 -a ../0_inputs/genomes/$speciesName.fa -o buscoAssemblyResults/$speciesName
            compleasm protein -l eukaryota -t 10 -p $species/protein.fasta -o buscoAnnotationResults/$speciesName
        fi
    done
done