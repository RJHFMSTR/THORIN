import gzip
import numpy as np
import time
from Bio import bgzf
from multiprocessing import Pool
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_vcf')
parser.add_argument('-o', '--output_prefix')
parser.add_argument('-p', '--pofo_probabilities')
args = parser.parse_args()

vcf_file=args.input_vcf
output_prefix=args.output_prefix
prob_file=args.pofo_probabilities






# Load probabilities into a dictionary
prob_dict = {}


with open(prob_file) as f:
    for line in f:
        if line.startswith('target'):
            continue
        tmp = line.split()
        prob_dict[tmp[0]]=float(tmp[4])










def parse_line(line: str, probs):

    initial_split = line.split('\t')
    fields = initial_split[:8]
    fields.append('GT:DS:GP')

    gt_fields = [x.split(":")[0] for x in initial_split[9:]]
    combined_ap_fields = [x.split(":")[2].partition(',') for x in initial_split[9:]]
    left_ap_fields = np.asarray([x[0] for x in combined_ap_fields], dtype=float)
    right_ap_fields = np.asarray([x[2] for x in combined_ap_fields], dtype=float)

    mat_ds = left_ap_fields * probs + right_ap_fields * (1 - probs)
    pat_ds = left_ap_fields * (1 - probs) + right_ap_fields * probs

    gt_fields_updated = ['.' if gt_fields[i] in ['0|0', '1|1'] else '0|0' if mat_ds[i] > pat_ds[i] else '1|1' for i in range(len(gt_fields))]

    diff_gp_str = ['./.:.:0,0,0' if gt_fields_updated[i] in ['.'] else f'{gt_fields_updated[i]}:{pat_ds[i]/(mat_ds[i]+pat_ds[i])}:{mat_ds[i]/(mat_ds[i]+pat_ds[i])},{pat_ds[i]/(mat_ds[i]+pat_ds[i])},0.0' for i in range(len(mat_ds))]
    diff_line = '\t'.join(fields + diff_gp_str)

    mat_gp_str = [f'{round(dsm)}|0:{dsm}:{1-dsm},{dsm},0' for dsm in mat_ds]
    mat_line = '\t'.join(fields + mat_gp_str)

    pat_gp_str = [f'{round(dsp)}|0:{dsp}:{1-dsp},{dsp},0' for dsp in pat_ds]
    pat_line = '\t'.join(fields + pat_gp_str)








    return [mat_line, pat_line, diff_line]









def process_vcf(vcf_file, PFX, prob_dict):


    with gzip.open(vcf_file, 'rt') as vcf_in, bgzf.open(PFX+'.maternal_haplotype.vcf.gz', 'wt') as out_mat, bgzf.open(PFX+'.paternal_haplotype.vcf.gz', 'wt') as out_pat, bgzf.open(PFX+'.differential_haplotype.vcf.gz', 'wt') as out_diff :
        for line in vcf_in:
            if line[0] == '#':
                out_mat.write(line)
                out_pat.write(line)
                out_diff.write(line)

                if line.find('#CHROM') == 0:
                    order_individual=line.split()[9:]
                    for t in order_individual:
                        if t not in prob_dict.keys():
                            prob_dict[t]=0.5
                    probs=np.asarray([prob_dict[t] for t in order_individual], dtype=float)

            else:
                out_lines  = parse_line(line, probs)

                out_mat.write(out_lines[0] + '\n')
                out_pat.write(out_lines[1] + '\n')
                out_diff.write(out_lines[2] + '\n')



start_time = time.time()
process_vcf(vcf_file, output_prefix, prob_dict)
end_time = time.time()
elapsed_time = end_time - start_time

print(f"Execution time: {elapsed_time:.2f} seconds")
