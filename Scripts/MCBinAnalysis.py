#!/usr/bin/env python3

import re
import csv
import math
import copy
import random

from argparse import ArgumentParser
from collections import defaultdict
from multiprocessing import Pool

import numpy as np

from Bio import SeqIO
from scipy.stats import norm, shapiro, false_discovery_control

def get_chromosome_length_and_bins(genome_fasta, bin_size):
    '''
    Use biopython to create a dictionary of chromosome lengths and lists for counts in each bin
    '''
    chrom_lengths = {}
    chrom_bins = {}

    for seq in SeqIO.parse(genome_fasta, 'fasta'):
        chrom_lengths[seq.id] = len(seq)
        chrom_bins[seq.id] = [0] * math.ceil(chrom_lengths[seq.id] / bin_size)

    return chrom_lengths, chrom_bins

def create_gene_bin_map(annotation_bed, bin_size):
    '''
    Create a map of gene ID -> [chromosome, bin] for all genes in annotation_bed
    '''
    gene_map = defaultdict(lambda: None)

    with open(annotation_bed) as f:
        reader = csv.reader(f, delimiter='\t')

        for line in reader:
            if line[7] == 'gene':
                gene_name = re.findall(r'ID=([^;]+)', line[-1])[0]
                gene_chrom = line[0]
                gene_centre = (int(line[1]) + int(line[2])) / 2
                gene_bin = int(gene_centre // bin_size)

                gene_map[gene_name] = [gene_chrom, gene_bin]

    return gene_map

def generate_random_distribution(n_genes):
    '''
    Select a random n_gene and calculate the distribution, for multiprocessing
    Accesses global bins and gene_to_bin_map
    '''
    dist_bins = copy.deepcopy(bins)

    for gene in random.sample(sorted(gene_to_bin_map.keys()), n_genes):
        if gene_to_bin_map[gene] is None:
            print(gene)
        chrom, bin_pos = gene_to_bin_map[gene]

        dist_bins[chrom][bin_pos] += 1

    return dist_bins

def count_dist_from_file(bins, gene_to_bin_map, gene_list_file):
    '''
    Read gene names from gene_list file and build a distribution of them
    '''
    dist = copy.deepcopy(bins)
    n_genes = 0

    with open(gene_list_file) as f:
        reader = csv.reader(f)

        for line in reader:
            if gene_to_bin_map[line[0]] is not None:
                gene_chrom, gene_bin = gene_to_bin_map[line[0]]
                dist[gene_chrom][gene_bin] += 1

                n_genes += 1
            else:
                del gene_to_bin_map[line[0]]

    return dist, n_genes

def calculate_ps_for_bins(bin_pos):
    '''
    Calculate the probablilty of a bin containing more genes than random chance
    Uses global random_dists and true_dist
    '''
    random_vals_for_bin = []
    true_value_for_bin = true_dist[bin_pos[0]][bin_pos[1]]
    for dist in random_dists:
        random_vals_for_bin.append(dist[bin_pos[0]][bin_pos[1]])

    mean, var = norm.fit(random_vals_for_bin)
    if true_value_for_bin < mean:
        return 1

    distrib = norm(mean, var)
    prob = distrib.pdf(true_value_for_bin)

    if np.isnan(prob):
        return 1
    else:
        return prob
    
def get_significance_marker(p_val):
    '''
    Return the significance marker for a p_val
    '''
    if p_val < 0.001:
        return '***'
    elif p_val < 0.01:
        return '**'
    elif p_val < 0.05:
        return '*'
    else:
        return ''

if __name__ == '__main__':
    parser = ArgumentParser(description='Run a monte carlo analysis on a gene set bin distribution and find bins with significantly more genes')

    parser.add_argument('genome_fasta', help='FASTA file containing the gneome assembly')
    parser.add_argument('annotation_bed', help='BED12 file containing the gene positions')
    parser.add_argument('selected_genes', help='File with one gene ID per line of interesting genes')

    parser.add_argument('-b', '--bin-size', help='Size in BP of bins', type=int, default=250_000)
    parser.add_argument('-r', '--random-runs', help='Number of random runs to perform to fit distributions', type=int, default=10_000)
    parser.add_argument('-o', '--output', help='File to write output probablities to', default='mc_results.tsv')
    parser.add_argument('-t', '--threads', help='Number of threads to use for multiprocessing', default=10, type=int)

    args = parser.parse_args()

    lengths, bins = get_chromosome_length_and_bins(args.genome_fasta, args.bin_size)
    gene_to_bin_map = create_gene_bin_map(args.annotation_bed, args.bin_size)

    true_dist, n_selected = count_dist_from_file(bins, gene_to_bin_map, args.selected_genes)

    distribution_inputs = [n_selected] * args.random_runs
    with Pool(args.threads) as p:
        random_dists = p.map(generate_random_distribution, distribution_inputs)

    bins_to_calc = []
    for chrom in bins:
        for i, value in enumerate(bins[chrom]):
            bins_to_calc.append([chrom, i])

    with Pool(args.threads) as p:
        p_vals = p.map(calculate_ps_for_bins, bins_to_calc)

    adj_ps = false_discovery_control(p_vals)

    output_file = [['Chromosome', 'Start', 'End', 'Bin Number', 'P value', 'Adj P', 'Significance']]
    for i, bin in enumerate(bins_to_calc):
        output_file.append([
            bin[0], bin[1] * args.bin_size, bin[1] * args.bin_size + args.bin_size - 1, bin[1], 
            p_vals[i], adj_ps[i], get_significance_marker(adj_ps[i])
        ])

    with open(args.output, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(output_file)