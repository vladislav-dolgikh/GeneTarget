import pandas as pd
import numpy
from Bio import SeqIO
from Bio import SeqFeature
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--chipseq_peaks', default='None', help = 'BED file contains ChIP-seq peaks')
parser.add_argument('--chipseq_fasta', help = 'FASTA file contains sequences of ChIP-seq peaks')
parser.add_argument('--rnaseq_table', help = 'Table contains log2FC and p-adj values from RNA-seq processing experiment')
parser.add_argument('--motif_recognition', help = 'A boolean function determines whether the motif will be recognized (0 - if not, 1 - if yes)')
parser.add_argument('--file_PWM', help = 'File contains position weight matrix')
parser.add_argument('--threshold', type = float, help = 'Threshold value for position weight matrix')
parser.add_argument('--file_output', help = 'Name of output file')
args = parser.parse_args()

def motif_finder(chipseq_peaks, file_PWM, chipseq_fasta, threshold):

    motif_recognition_need = chipseq_peaks

    if str(motif_recognition_need) == 'None':

        file_temp = open('temp.txt', 'a')

        for seq_record in SeqIO.parse(open(chipseq_fasta), "fasta"):
            peak_coordinate = seq_record.id
            peak_coordinate = peak_coordinate.split('-')
            peak_chromosome = str(peak_coordinate[0])
            peak_start = str(peak_coordinate[1])
            peak_end = str(peak_coordinate[2])
            file_temp.write(str(peak_chromosome) + '\t' + str(peak_start) + '\t' + str(peak_end) + '\n')

    else:

        # Upload PWM file, entering the values of motif length and motif threshold
        file_temp = open('temp.txt', 'a')
        file_PWM=numpy.loadtxt(open(file_PWM, 'r'))
        motif_length = int(len(file_PWM))
        file_PWM=numpy.transpose(file_PWM)
        motif_threshold=float(threshold)

        # Counting of minimum and maximum weight of PWM
        wmax = 0
        wmin = 0
        for i in range (0, motif_length):
            wmax = wmax + max(file_PWM[:,i])
            wmin = wmin + min(file_PWM[:,i])

        # making the reverse complement strand of fasta file
        backup = open("backup.txt","a")
        for seq_record in SeqIO.parse(open(chipseq_fasta), "fasta"):
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
                    peak_coordinate = peak_coordinate.split('-')
                    peak_chromosome = str(peak_coordinate[0])
                    peak_start = str(peak_coordinate[1])
                    peak_end = str(peak_coordinate[2])
                    file_temp.write(str(peak_chromosome) + '\t' + str(peak_start) + '\t' + str(peak_end) + '\t' + str(feature) + '\n')
                g += 1
                t += 1
                v += 1
        backup.close()
        file_temp.close()

def peak_in_promoters(chipseq_peaks):

    chipseq_peaks = open(chipseq_peaks, "r")  # BED file with ChIP-seq peaks
    for line1 in chipseq_peaks:
        line1 = line1.strip('\n')
        line1 = line1.split('\t')
        peak_start = int(line1[1])
        peak_end = int(line1[2])
        chipseq_chromosome = line1[0]
        if len(line1) == 4:
            site = str(line1[3])
        promoters_coordinates = open("volleydoll_promoters_tair10.txt", "r")  # BED file with promoters coordinates
        for line2 in promoters_coordinates:
            line2 = line2.strip('\n')
            line2 = line2.split('\t')
            region_start = int(line2[1])
            region_end = int(line2[2])
            promoters_chromosome = line2[0]
            gene_name = line2[3]
            if chipseq_chromosome == promoters_chromosome:
                if max(peak_start, region_start) < min(peak_end, region_end):
                    file3 = open('test_intersect.txt', "a")  # output
                    start = max(peak_start, region_start)
                    end = min(peak_end, region_end)
                    length = min(peak_end, region_end) - max(peak_start, region_start)
                    if len(line1) == 3:
                        file3.write(str(chipseq_chromosome) + '\t' + str(start) + '\t' + str(end) + '\t' + str(gene_name) + '\n')
                    elif len(line1) == 4:
                        file3.write(str(chipseq_chromosome) + '\t' + str(start) + '\t' + str(end) + '\t' + str(gene_name) + '\t' +
                                    str(site) + '\n')
                    file3.close()
        promoters_coordinates.close()
    chipseq_peaks.close()

def target_gene(rnaseq_table, file_output):
    rna_seq = pd.read_csv(rnaseq_table, sep="\t", header=None)
    targets_genelist = pd.read_csv('test_intersect.txt', sep="\t", header=None)

    experiments_count = int((rna_seq.shape[1] - 1) / 2)
    gene_number = int(rna_seq.shape[0])
    target_number = int(targets_genelist.shape[0])
    site_presence = int(targets_genelist.shape[1]) # 4 - no; 5 - yes

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
                a.append('non_differentially_expressed') # С‡РµС‚РЅС‹Р№ СЃС‚РѕР»Р±РµС†
            j = j + 1


        if 'differentially_expressed' in a and 'up-regulated' in a:
            up_regulated_degs.append(gene_attributes.tolist())
        elif 'differentially_expressed' in a and 'down-regulated' in a:
            down_regulated_degs.append(gene_attributes.tolist())

        i = i + 1

    output = open(file_output, "a")

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
                    if site_presence == 4:
                        output.write(str(targets_genelist.iat[e,0]) + '\t' + str(targets_genelist.iat[e,1]) + '\t' + str(targets_genelist.iat[e,2]) + '\t' +
                                     str(targets_genelist.iat[e, 3]) + '\t')
                        output.writelines("%s\t" % line for line in deg_attribute)
                        output.write('\n')
                    elif site_presence == 5:
                        output.write(str(targets_genelist.iat[e, 0]) + '\t' + str(targets_genelist.iat[e, 1]) + '\t' + str(targets_genelist.iat[e, 2]) + '\t' +
                            str(targets_genelist.iat[e, 3]) + '\t' + str(targets_genelist.iat[e, 4]) + '\t')
                        output.writelines("%s\t" % line for line in deg_attribute)
                        output.write('\n')
                e = e + 1
            d = d + 1
    annotator(up_regulated_degs)
    annotator(down_regulated_degs)

    uniqlines = set(open(file_output, 'r', encoding='utf-8').readlines())
    final = open(file_output, 'w', encoding='utf-8').writelines(set(uniqlines))


motif_finder(chipseq_peaks = args.chipseq_peaks, file_PWM = args.file_PWM, chipseq_fasta = args.chipseq_fasta, threshold = args.threshold)
peak_in_promoters(chipseq_peaks = 'temp.txt')
target_gene(rnaseq_table = args.rnaseq_table, file_output = args.file_output)

os.remove('backup.txt')
os.remove('temp.txt')
os.remove('test_intersect.txt')
