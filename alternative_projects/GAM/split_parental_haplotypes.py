import argparse
import gzip
import os
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_vcf')
parser.add_argument('-h1', '--output_haplotype1')
parser.add_argument('-h2', '--output_haplotype2')
args = parser.parse_args()

vcf_in=args.input_bcf
out_h1=args.output_haplotype1
out_h2=args.output_haplotype2

outfile1=open(out_h1,'w')
outfile2=open(out_h2,'w')
with gzip.open(vcf_in,'rb') as f:
    for line1 in f:
        line=line1.decode()
        tmp=line.split()
        if line.find('#')==0:
            outfile1.write(line)
            outfile2.write(line)

        else:
            w1=tmp[:8]; w1.append('GT:DS:GP')
            w2=tmp[:8]; w2.append('GT:DS:GP')
            for fields in tmp[9:]:

                GT=fields.split(':')[0]
                GT1=GT.split('|')[0]
                GT2=GT.split('|')[1]

                AP=fields.split(':')[2]
                AP1=AP.split(',')[0]
                AP2=AP.split(',')[1]

                GP1=str((1-float(AP1)))+','+AP1+',0'
                GP2=str((1-float(AP2)))+','+AP2+',0'

                field_h1=GT1+'|0:'+AP1+':'+GP1
                field_h2=GT2+'|0:'+AP2+':'+GP2    

                w1.append(field_h1)
                w2.append(field_h2)

            outfile1.write('\t'.join(w1)+'\n')
            outfile2.write('\t'.join(w2)+'\n')

outfile1.close()
outfile2.close()


os.system('bgzip -f '+out_h1)
os.system('bgzip -f '+out_h2)
os.system('bcftools index -f '+out_h1+'.gz')
os.system('bcftools index -f '+out_h2+'.gz')
