#!/usr/bin/env python3

import os
import csv
import glob

from argparse import ArgumentParser
from collections import defaultdict
from statistics import mean

from Bio import SeqIO
from Bio.SeqUtils import GC

from pronto import Ontology
from gffutils import DataIterator
from scipy.stats import mannwhitneyu

SCRIPT_PATH = os.path.realpath(os.path.dirname(__file__))
GO_OBO = os.path.join(SCRIPT_PATH, 'nameMapFiles', 'go.obo')
KEGG_MAP = os.path.join(SCRIPT_PATH, 'nameMapFiles', 'KEGGNameMap.tsv')
PFAM_MAP = os.path.join(SCRIPT_PATH, 'nameMapFiles', 'PfamNameMap.tsv')

def get_significance_mark(p_value):
    '''
    Return stars equal to significance of a p value
    '''
    if p_value < 0.001:
        return '***'
    elif p_value < 0.01:
        return '**'
    elif p_value < 0.05:
        return '*'
    else:
        return ''

def get_genome_stats(genome_fasta):
    '''
    calculate genome length, N50 and number of scaffolds
    '''
    scaffold_lengths = []
    num_scaffolds = 0
    gc_bases = 0

    for seq in SeqIO.parse(genome_fasta, 'fasta'):
        num_scaffolds += 1
        scaffold_lengths.append(len(seq))

        gc_bases += len(seq) * GC(seq.seq)

    scaffold_lengths.sort(reverse=True)
    total_length = sum(scaffold_lengths)
    half_length = total_length / 2

    size_total = 0
    n50 = 0
    for scaff_size in scaffold_lengths:
        n50 += 1
        size_total += scaff_size

        if size_total > half_length:
            break

    return total_length, n50, gc_bases / total_length, num_scaffolds

    
def parse_island_info(islands_bed, islands_prob):
    '''
    Parse the number and expected number of islands from a bed and txt file
    '''
    num_islands = 0
    num_genes_in_islands = 0
    genes_in_islands = defaultdict(lambda: [])
    island_lengths = defaultdict(lambda: 0)

    with open(islands_bed) as f:
        reader = csv.reader(f, delimiter='\t')

        for line in reader:
            num_islands += 1
            genes = line[-1].split(',')

            num_genes_in_islands += len(genes)
            genes_in_islands[f'{line[0]}:{line[1]}-{line[2]}'] = genes
            island_lengths[f'{line[0]}:{line[1]}-{line[2]}'] = int(line[2]) - int(line[1])

    prob = 1
    significance = ''
    random_num_clusters = []
    expected_genes = []
    expected_lengths = []
    with open(islands_prob) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)

        for line in reader:
            if len(line) > 0 and line[0] == 'Probablity of number of clusters (Generalised Poisson):':
                prob = float(line[1])
                significance = line[2]

                break
            elif len(line) == 5:
                random_num_clusters.append(int(line[0]))
                genes = line[4].split(',')
                lengths = line[2].split(',')

                if genes[0] != '':
                    expected_genes.append(list(map(int, genes)))

                if lengths[0] != '':
                    expected_lengths.append(list(map(int, lengths)))

    if len(random_num_clusters) == 0:
        expected_clusters = 0
    else:
        expected_clusters = sum(random_num_clusters) / len(random_num_clusters)

    return num_islands, expected_clusters, prob, significance, genes_in_islands, island_lengths, expected_genes, expected_lengths

def read_busco_score(busco_output_dir):
    '''
    Read the BUSCO completeness from compleasm output
    '''
    results = {}
    with open(os.path.join(busco_output_dir, 'summary.txt')) as f:
        reader = csv.reader(f, delimiter=':')
        next(reader)

        for line in reader:
            results[line[0]] = float(line[1].split('%')[0])

    return results['S'] + results['D']

def read_name_clade_host(tsv_file):
    '''
    Read information on the name clade and host of the species from the TSV file
    '''
    species = []
    with open(tsv_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)

        for line in reader:
            species.append(line)

    return species

def read_name_map(name_map_file):
    '''
    Read a name map into a defaultdict
    '''
    name_map = defaultdict(lambda: 'NONE')

    with open(name_map_file) as f:
        reader = csv.reader(f, delimiter='\t')

        for line in reader:
            name_map[line[0]] = line[1]

    return name_map

def read_enriched_fams(enrichment_file, go_name_map, kegg_name_map, pfam_name_map):
    '''
    Read the enriched GO terms, kegg pathways and PFams from the TSVs
    '''
    enrichment_lines = {}
    organism_enrichment = []
    with open(enrichment_file) as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader)

        for line in reader:
            new_line = [line[1], line[2], line[3]]

            # Add names to GO Terms
            go_names = []
            for go_term in line[4].split(','):
                if go_term != '':
                    try:
                        go_names.append(go_name_map[go_term].name)
                    except:
                        go_names.append('NONE')

            new_line.append(line[4])
            new_line.append(','.join(go_names))
            new_line.append(line[5])

            kegg_names = []
            for kegg in line[6].split(','):
                if kegg != '':
                    kegg_names.append(kegg_name_map[kegg])

            new_line.append(line[6])
            new_line.append(','.join(kegg_names))

            kegg_names = []
            for kegg in line[7].split(','):
                if kegg != '':
                    kegg_names.append(kegg_name_map[kegg])

            new_line.append(line[7])
            new_line.append(','.join(kegg_names))
            new_line.append(line[8])
            new_line.append(line[9])

            pfam_names = []
            for pfam in line[10].split(','):
                if pfam != '':
                    pfam_names.append(pfam_name_map[pfam])

            new_line.append(line[10])
            new_line.append(','.join(pfam_names))

            pfam_names = []
            for pfam in line[11].split(','):
                if pfam != '':
                    pfam_names.append(pfam_name_map[pfam])

            new_line.append(line[11])
            new_line.append(','.join(pfam_names))
            new_line.append(line[12])
            new_line.append(line[13])

            if line[0] in ['all_islands', 'all_inverse']:
                organism_enrichment = new_line
            else:
                enrichment_lines[line[0]] = new_line

    return organism_enrichment, enrichment_lines

def count_genes(gff_file):
    '''
    Count the number of genes for a species from it's GFF file
    '''
    gene_count = 0
    for feature in DataIterator(gff_file):
        if feature.featuretype == 'gene':
            gene_count += 1

    return gene_count

def calc_length_num_significance(true_vals, expected_vals):
    '''
    Calculate averages and significance of true vals vs expected
    '''
    if len(expected_vals) > 0 and type(expected_vals[0]) == list:
        expected_vals = [x for xs in expected_vals for x in xs]

    if len(true_vals) > 0 and type(true_vals[0]) == list:
        true_vals = [x for xs in true_vals for x in xs]

    if len(true_vals) > 0 and len(expected_vals) > 0: 
        true_av = mean(true_vals)
        expected_av = mean(expected_vals)

        mwu = mannwhitneyu(true_vals, expected_vals)
        p_val = mwu.pvalue

        if p_val >= 0.05:
            differance = 'Not Sig'
        elif true_av > expected_av:
            differance = 'Greater'
        else:
            differance = 'Less'

        return true_av, expected_av, differance, p_val
    elif len(true_vals) > 0:
        expected_av = mean(expected_vals)

        return 0, expected_av, 'Not Sig', 1
    elif len(expected_vals) > 0:
        true_av = mean(true_vals)

        return true_av, 0, 'Not Sig', 1
    else:
        return 0, 0, 'Not Sig', 1

if __name__ == '__main__':
    parser = ArgumentParser(description='Create supplementry tables for the species overview and parasitism islands')

    parser.add_argument('name_clade_list', help='Path to the TSV containing name clade and host info')
    parser.add_argument('genomes_directory', help='Directory containing the genomes for the analysis')
    parser.add_argument('results_directory', help='The directory containing the pipeline results')
    parser.add_argument('busco_assembly_out', help='Directory containing the compleasm BUSCO assembly output')
    parser.add_argument('busco_annotation_out', help='Directory containing the compleasm BUSCO annotation output')
    parser.add_argument('enriched_term_dir', help='Path to the directory with the enriched terms')

    parser.add_argument('-o', '--output', help='TSV file to write output table to', default='suppTable')
    parser.add_argument('-t', '--threads', help='Number of threads to use to run this', default=10, type=int)

    args = parser.parse_args()

    summary_table = [[
        'Name', 'Clade', 'Host', 'Genome Number of Scaffolds', 'Genome N50', 'Genome Size (bp)', 'BUSCO Assembly', 'Genome-Wide GC content (%)',
        'Total Number of Annotated Genes', 'BUSCO Annotation', 'Number of Islands', 'Expected Number of Islands', 'Number of Islands P Value', 'Number of Islands P Value Significance',
        'Mean Number of Genes in Islands', 'Expected Number of Genes in Islands', 'Direction', 'Number of Genes P Value', 'Number of Genes P Value Significance',
        'Mean Length (bp) of Island', 'Expected Length (bp) of Islands', 'Direction', 'Length of Islands P Value', 'Length of Islands P Value Significance',
        'Number of Inverse Islands', 'Expected Number of Inverse Islands', 'Number of Inverse Islands P Value', 'Number of Inverse Islands P Value Significance'
    ]]

    species_enrichment_table = [[
        'Species', 'eggNOG Group Count', 'All GO IDs In Islands', 'Enriched GO IDs in Island', 'Enriched GO Terms in Island',
        'TopGO P Values', 'All KEGG IDs in Island', 'All KEGG Pathways in Island', 'Enriched KEGG IDs in Island', 'Enriched KEGG Pathways in Island', 
        'KEGG P Values', 'KEGG FDRs', 'All PFams in Island (Short Names)', 'All PFam in Island (Long Names)', 'Enriched PFam in Island (Short Names)', 
        'Enriched PFam in Island (Long Names)', 'PFam P Values', 'PFam FDRs', 'eggNOG Group Count', 'All GO IDs In Islands',
        'Enriched GO IDs in Island', 'Enriched GO Terms in Island',
        'TopGO P Values', 'All KEGG IDs in Island', 'All KEGG Pathways in Island', 'Enriched KEGG IDs in Island', 'Enriched KEGG Pathways in Island', 
        'KEGG P Values', 'KEGG FDRs', 'All PFams in Island (Short Names)', 'All PFam in Island (Long Names)', 'Enriched PFam in Island (Short Names)', 
        'Enriched PFam in Island (Long Names)', 'PFam P Values', 'PFam FDRs'
    ]]

    island_info_table = [[
        'Species', 'Island Number', 'Chromosome', 'Start', 'End', 'Length (bp)', 'Number of Genes', 'Gene IDs',
        'eggNOG Group Count', 'Island Type', 'Enriched GO IDs in Island', 'Enriched GO Terms in Island',
        'TopGO P Values', 'All KEGG IDs in Island', 'All KEGG Pathways in Island', 'Enriched KEGG IDs in Island', 'Enriched KEGG Pathways in Island', 
        'KEGG P Values', 'KEGG FDRs', 'All PFams in Island (Short Names)', 'All PFam in Island (Long Names)', 'Enriched PFam in Island (Short Names)', 
        'Enriched PFam in Island (Long Names)', 'PFam P Values', 'PFam FDRs'
    ]]

    inverse_island_info_table = [[
        'Species', 'Inverse Island Number', 'Chromosome', 'Start', 'End', 'Length (bp)', 'Number of Genes', 'Gene IDs',
        'eggNOG Group Count', 'Inverse Island Type', 'All GO IDs In Islands', 'Enriched GO IDs in Island', 'Enriched GO Terms in Island',
        'TopGO P Values', 'All KEGG IDs in Island', 'All KEGG Pathways in Island', 'Enriched KEGG IDs in Island', 'Enriched KEGG Pathways in Island', 
        'KEGG P Values', 'KEGG FDRs', 'All PFams in Island (Short Names)', 'All PFam in Island (Long Names)', 'Enriched PFam in Island (Short Names)', 
        'Enriched PFam in Island (Long Names)', 'PFam P Values', 'PFam FDRs'
    ]]

    go_obo = Ontology(GO_OBO)
    kegg_name_map = read_name_map(KEGG_MAP)
    pfam_name_map = read_name_map(PFAM_MAP)

    for species, clade, host in read_name_clade_host(args.name_clade_list):
        print(f'==> Extracting data for {species}...')

        species_underscores = species.replace(' ', '_')
        species_dir = os.path.dirname(glob.glob(os.path.join(args.results_directory, '*', species_underscores, 'clusters.dbscan.bed'))[0])
        genome_fasta = os.path.join(args.genomes_directory, f'{species_underscores}.fa')

        genome_stats = get_genome_stats(genome_fasta)
        island_info = parse_island_info(
            os.path.join(species_dir, 'clusters.dbscan.bed'),
            os.path.join(species_dir, 'clusterStatistics.dbscan.tsv')
        )
        inverse_island_info = parse_island_info(
            os.path.join(species_dir, 'clusters_inverse.dbscan.bed'),
            os.path.join(species_dir, 'clusterStatisticsDownreg.dbscan.tsv')
        )
        
        length_stats = calc_length_num_significance(list(island_info[5].values()), island_info[7])

        gene_nums = []
        for cluster in island_info[4]:
            gene_nums.append(len(island_info[4][cluster]))
        num_stats = calc_length_num_significance(gene_nums, island_info[6])

        num_genes = count_genes(os.path.join(species_dir, 'filtered.gff3'))

        busco_assembly = read_busco_score(os.path.join(
            args.busco_assembly_out, species_underscores
        ))
        busco_annotation = read_busco_score(os.path.join(
            args.busco_annotation_out, species_underscores
        ))

        # Name	Clade	Lifestyle	Genome Number of Scaffolds	Genome N50
        # Genome Size (bp)	BUSCO Assembly	Genome-Wide GC content (%)
        # Total Number of Annotated Genes	BUSCO Annotation	Average Number of Genes in Islands
        # Expected Number of Genes in Islands	Average Length of Islands	Expected Length of Islands
        # Average Island GC Content	Number of Clusters	Expected Number of Clusters	P Value
        # P Value Significance	Number of Inverse Clusters	Expected Number of Inverse Clusters	P Value
        # P Value Significance
        summary_table.append([
            species, clade, host, genome_stats[3], genome_stats[1], genome_stats[0], round(busco_assembly, 2), genome_stats[2],
            num_genes, busco_annotation, 
            island_info[0], island_info[1], island_info[2], get_significance_mark(island_info[2]),
            num_stats[0], num_stats[1], num_stats[2], num_stats[3], get_significance_mark(num_stats[3]),
            length_stats[0], length_stats[1], length_stats[2], length_stats[3], get_significance_mark(length_stats[3]),
            inverse_island_info[0], inverse_island_info[1], inverse_island_info[2], get_significance_mark(inverse_island_info[2])
        ])

        island_enrichment = read_enriched_fams(os.path.join(args.enriched_term_dir, f'{species_underscores}_islands.tsv'), go_obo, kegg_name_map, pfam_name_map)
        inverse_island_enrichment = read_enriched_fams(os.path.join(args.enriched_term_dir, f'{species_underscores}_inverse.tsv'), go_obo, kegg_name_map, pfam_name_map)

        # Species	Cluster Number	Chromosome	Start	End	Length (bp)	Number of Genes	Gene ID List
        # EggNOG Group Count  Island Type Classification
        # GO IDs	GO Terms	GO P Values	KEGG IDs	KEGG Names	KEGG P Values	PFams
        # PFam P Values	Interpro Domains	Interpro Names	Interpro P Values	
        for i, island in enumerate(island_info[4]):
            position = island.split(':')
            chrom = position[0]
            coords = position[1].split('-')
            start = int(coords[0])
            end = int(coords[1])

            island_info_table.append([
                species, i + 1, chrom, start, end, end - start, len(island_info[4][island]), ','.join(island_info[4][island])
            ] + island_enrichment[1][island])

        for i, island in enumerate(inverse_island_info[4]):
            position = island.split(':')
            chrom = position[0]
            coords = position[1].split('-')
            start = int(coords[0])
            end = int(coords[1])

            inverse_island_info_table.append([
                species, i + 1, chrom, start, end, end - start, len(inverse_island_info[4][island]), ','.join(inverse_island_info[4][island])
            ] + inverse_island_enrichment[1][island])

        species_enrichment_table.append([species] + [island_enrichment[0][0]] + island_enrichment[0][2:] + [inverse_island_enrichment[0][0]] + inverse_island_enrichment[0][2:])

    print('==> Writing tables...')

    with open('summaryTable.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(summary_table)

    with open('islandTable.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(island_info_table)

    with open('inverseIslandTable.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(inverse_island_info_table)

    with open('speciesEnrichment.tsv', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerows(species_enrichment_table)