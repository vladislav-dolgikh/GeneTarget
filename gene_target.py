import pandas as pd
import numpy
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('chipseq_peaks', help = 'BED file contains ChIP-seq peaks')
parser.add_argument('rnaseq_table', help = 'Table contains log2FC and p-adj values from RNA-seq processing experiment')
parser.add_argument('motif_recognition', help = 'A boolean function determines whether the motif will be recognized (0 - if not, 1 - if yes)')
parser.add_argument('file_PWM', help = 'File contains position weight matrix')
parser.add_argument('threshold', type = float, help = 'Threshold value for position weight matrix')
parser.add_argument('chipseq_fasta', help = 'FASTA file contains sequences of ChIP-seq peaks')
parser.add_argument('file_output', help = 'Name of output file')
args = parser.parse_args()

def motif_finder():

    # Upload PWM file, entering the values of motif length and motif threshold
    file_temp = open('temp.txt', 'a')
    file_PWM=numpy.loadtxt(open(args.file_PWM, 'r'))
    motif_length = int(len(file_PWM))
    file_PWM=numpy.transpose(file_PWM)
    motif_threshold=float(args.threshold)

    # Counting of minimum and maximum weight of PWM
    wmax = 0
    wmin = 0
    for i in range (0, motif_length):
        wmax = wmax + max(file_PWM[:,i])
        wmin = wmin + min(file_PWM[:,i])

    # making the reverse complement strand of fasta file
    backup = open("backup.txt","a")
    for seq_record in SeqIO.parse(open(args.promoter_sequence), "fasta"):
        peak = seq_record.seq
        peak_rev = peak.reverse_complement()
        backup.write(str(seq_record.id) + '\t' + str(peak) + '\t' + str('forward') + '\n' +
                     str(seq_record.id) + '\t' + str(peak_rev) +  '\t' + str('reverse') + '\n')
    backup.close()

    # motif finding
    peaks = open("backup.txt","r")
    for line in peaks:
        line = line.strip('\n')
        line = line.split('\t')
        id = line[0]
        sequence = line[1]
        orientation = line[2]
        len_sequence = len(sequence)
        g = 0
        t = motif_length
        v = 0
        while v < (len_sequence - (motif_length - 1)):
            peak_coordinate = id
            feature = sequence[g:t]
            score = 0
            for i in range(0, motif_length):
                if feature[i] == 'A':
                    score += (file_PWM[0, i])
                elif feature[i] == 'C':
                    score += (file_PWM[1, i])
                elif feature[i] == 'G':
                    score += (file_PWM[2, i])
                elif feature[i] == 'T':
                    score += (file_PWM[3, i])
                elif feature[i] == 'n':
                    score += 100000
            weighted_score = (score - wmin) / (wmax - wmin)
            if weighted_score >= motif_threshold and weighted_score <= 1:
                if orientation == 'forward':
                    peak_coordinate = peak_coordinate.split('-')
                    motif_coordinate_chromosome = str(peak_coordinate[0])
                    motif_coordinate_start = (int(peak_coordinate[1])) + g
                    motif_coordinate_end = (int(peak_coordinate[1])) + t
                    file_temp.write(str(peak_coordinate) + '\t' + str(feature) + '\t' + str(weighted_score) + '\t' + str(motif_coordinate_chromosome) + '\t'
                                      + str(motif_coordinate_start) + '\t' + str(motif_coordinate_end) + '\t' + str('+') + '\n')
                if orientation == 'reverse':
                    peak_coordinate = peak_coordinate.split('-')
                    motif_coordinate_chromosome = str(peak_coordinate[0])
                    motif_coordinate_end = (int(peak_coordinate[2])) - g
                    motif_coordinate_start = (int(peak_coordinate[2])) - t
                    file_temp.write(str(peak_coordinate) + '\t' + str(feature) + '\t' + str(weighted_score) + '\t' + str(motif_coordinate_chromosome) + '\t'
                                      + str(motif_coordinate_start) + '\t' + str(motif_coordinate_end) + '\t' + str('-') + '\n')
            g += 1
            t += 1
            v += 1
    backup.close()
    file_temp.close()

    promoters_coordinates = open(args.promoter_coordinates, 'r')
    for line1 in promoters_coordinates:
        line1 = line1.strip('\n')
        line1 = line1.split('\t')
        c1 = int(line1[0])
        x1 = int(line1[1])
        y1 = int(line1[2])
        gene_name = line1[3]
        promoter_strand = line1[4]
        found_sites = open('temp.txt', 'r')
        for line2 in found_sites:
            line2 = line2.strip('\n')
            line2 = line2.split('\t')
            c2 = int(line2[3])
            x2 = int(line2[4])
            y2 = int(line2[5])
            site = line2[1]
            peak = line2[0]
            site_strand = line2[6]
            if c1 == c2:
                if max(x1, x2) < min(y1, y2) and promoter_strand == '+':
                    file_out = open('temp_2.txt', 'a')
                    start = max(x1, x2)
                    end = min(y1, y2)
                    length = min(y1, y2) - max(x1, x2)
                    to_tss = y1 - end
                    print(start)
                    print(x2)
                    print(to_tss)
                    file_out.write(str(c1) + '\t' + str(start) + '\t' + str(end) + '\t' + str(to_tss) + '\t'
                                   + str(length) + '\t' + str(site) + '\t' + str(gene_name) + '\t'
                                   + str(site_strand) + '\t' + str(peak) + '\n')
                elif max(x1, x2) < min(y1, y2) and promoter_strand == '-':
                    file_out = open('temp_2.txt', 'a')
                    start = max(x1, x2)
                    end = min(y1, y2)
                    length = min(y1, y2) - max(x1, x2)
                    to_tss = end - x1
                    print(start)
                    print(x2)
                    print(to_tss)
                    file_out.write(str(c1) + '\t' + str(start) + '\t' + str(end) + '\t' + str(to_tss) + '\t'
                                   + str(length) + '\t' + str(site) + '\t' + str(gene_name) + '\t'
                                   + str(site_strand) + '\t' + str(peak) + '\n')

    found_sites.close()
    promoters_coordinates.close()
    file_out.close()

    annotated_genes = pd.read_csv('temp_2.txt', sep="\t", header=0,
                                  names=['Chromosome', 'Start', 'End', 'ToTSS', 'Length', 'Site', 'Gene', 'Orientation',
                                         'Peak'])
    annotated_genes = annotated_genes.sort_values(axis=0, by=['Chromosome', 'Start', 'End', 'ToTSS'],
                                                  ascending=[True, True, True, True])
    annotated_genes.to_csv('temp_3.txt', sep='\t', header=False, index=False)

    input_1 = pd.read_csv('temp_3.txt', sep="\t", header=None)

    tags = []
    tags.append(int(3000))
    i = 0
    while i < int(len(input_1) - 1):
        a = input_1.iat[i, 1]
        b = input_1.iat[i + 1, 1]
        tags.append(int(b - a))
        i = i + 1
    input_1['9'] = tags
    input_1.to_csv('temp_4.txt', sep='\t', header=False, index=False)

    input_2 = open('temp_4.txt', 'r')
    output_1 = open('temp_5.txt', 'a')
    for line in input_2:
        line = line.strip('\n')
        line = line.split('\t')
        tag = line[9]
        if int(tag) > 12 or int(tag) < 0:
            output_1.write(
                str(line[0]) + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\t' + str(line[4])
                + '\t' + str(line[5]) + '\t' + str(line[6]) + '\t' + str(line[7]) + '\t' + str(line[8]) + '\n')

def peak_in_promoters():

    chipseq_peaks = open(args.chipseq_peaks, "r")  # BED file with ChIP-seq peaks

    for line1 in chipseq_peaks:
        line1 = line1.strip('\n')
        line1 = line1.split('\t')
        peak_start = int(line1[1])
        peak_end = int(line1[2])
        chipseq_chromosome = line1[0]

        promoters_coordinates = open("promoters_tair10.txt", "r")  # BED file with promoters coordinates
        for line2 in promoters_coordinates:
            line2 = line2.strip('\n')
            line2 = line2.split('\t')
            region_start = int(line2[1])
            region_end = int(line2[2])
            promoters_chromosome = line2[0]
            gene_name = line2[3]
            if chipseq_chromosome == promoters_chromosome:
                if max(peak_start, region_start) < min(peak_end, region_end):
                    file3 = open('test_intersect_3.txt', "a")  # output
                    start = max(peak_start, region_start)
                    end = min(peak_end, region_end)
                    length = min(peak_end, region_end) - max(peak_start, region_start)
                    file3.write(str(chipseq_chromosome) + '\t' + str(start) + '\t' + str(end) + '\t' + str(gene_name) + '\n')
                    file3.close()
        promoters_coordinates.close()
    chipseq_peaks.close()

def target_gene():
    rna_seq = pd.read_csv('rnaseq_test.txt', sep="\t", header=None)
    targets_genelist = pd.read_csv('test_intersect_3.txt', sep="\t", header=None)

    experiments_count = int((rna_seq.shape[1] - 1) / 2)
    gene_number = int(rna_seq.shape[0])
    target_number = int(targets_genelist.shape[0])

    i = 0
    up_regulated_degs = list()
    down_regulated_degs = list()
    while i < gene_number:
        j = 1
        gene_attributes = rna_seq.iloc[i]
        lst = gene_attributes.shape[0]
        a = list()
        a.append(gene_attributes[j - 1])
        while j < lst:
            if j%2 != 0 and gene_attributes[j] > 0:
                a.append('up-regulated')
            elif j%2 != 0 and gene_attributes[j] < 0:
                a.append('down-regulated')
            elif j%2 == 0 and gene_attributes[j] < 0.05:
                a.append('differentially_expressed')
            elif j%2 == 0 and gene_attributes[j] > 0.05:
                a.append('non_differentially_expressed') # четный столбец
            j = j + 1


        if 'differentially_expressed' in a and 'up-regulated' in a:
            print('yes')
            up_regulated_degs.append(gene_attributes.tolist())
        elif 'differentially_expressed' in a and 'down-regulated' in a:
            print('no')
            down_regulated_degs.append(gene_attributes.tolist())
        print(up_regulated_degs)
        print('\n')
        print(down_regulated_degs)

        i = i + 1

    output = open("output.bed", "a")

    def annotator(deg_list):
        deg_list = deg_list
        d = 0
        while d < len(deg_list):
            e = 0
            while e < len(targets_genelist):
                deg_attribute = deg_list[d]
                deg_gene_name = deg_attribute[0]
                target_gene_name = targets_genelist.iat[e,3]
                if deg_gene_name == target_gene_name:
                    output.write(str(targets_genelist.iat[e,0]) + '\t' + str(targets_genelist.iat[e,1]) + '\t' + str(targets_genelist.iat[e,2]) + '\t' +
                                 str(targets_genelist.iat[e, 3]) + '\t')
                    output.writelines("%s\t" % line for line in deg_attribute)
                    output.write('\n')
                print('\n')
                e = e + 1
            d = d + 1
    annotator(up_regulated_degs)
    annotator(down_regulated_degs)

    uniqlines = set(open("output.bed", 'r', encoding='utf-8').readlines())
    final = open("output.bed", 'w', encoding='utf-8').writelines(set(uniqlines))

def launch():
    if args.motif_recognition == 0:
        peak_in_promoters()
        target_gene()
    elif args.motif_recognition == 1:
        motif_finder()
        peak_in_promoters()
        target_gene()

launch()
