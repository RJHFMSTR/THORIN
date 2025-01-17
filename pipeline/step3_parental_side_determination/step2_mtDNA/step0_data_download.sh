#!bin/bash

# Input file containing the CRAM URLs
index_url=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/1000genomes.high_coverage.GRCh38DH.alignment.index
index_file=$(basename ${index_url})
wget ${index_url}

# Output directory for chrM BAM files
ODIR=data/chrM_bams
mkdir -p ${ODIR}

# Extract CRAM URLs (ignoring header lines)
awk 'NR > 8 {print $1}' ${index_file} | sed 's|ftp:/|ftp://|' > cram_urls.txt

# Reference data
REF=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

# Download and extract chrM for each CRAM file
while read -r cram_url; do
    # Extract sample name from CRAM URL
    sample_name=$(basename "$cram_url" .cram)

    # Path to the output BAM file
    output_bam=${ODIR}/${sample_name}_chrM.bam

    # Download and extract chrM using samtools
    echo "Processing ${sample_name}..."
    samtools view -b -T ${REF} ${cram_url} chrM > "$output_bam"

    if [[ $? -eq 0 ]]; then
        echo "chrM extracted for ${sample_name} and saved to ${output_bam}"
    else
        echo "Failed to process ${sample_name}"
    fi
done < cram_urls.txt

rm cram_urls.txt
rm ${index_file}

