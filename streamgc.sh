#!/bin/bash

# Print help message if --help is called
if [[ "$1" == "--help" ]]; then
    echo "Program: streamGC"
    echo "Version: 1.0.0"
    echo "Code:    https://github.com/yh1126611/streamGC"
    echo "Usage:   streamgc <genome.fasta> <coordinates.txt> <output.tsv>"
    echo ""
    echo "Calculate the GC content (0-1) for every 100 bp window inside 10,000 bp interval from set of coordinates on genome."
    echo "<coordinates.txt>: a three-column tab-delimited file. Column 1: chromosome, column 2: genomic coordinate, column 3: strand orientation (+|-)"
    exit 0
fi

# Check if required arguments are provided
if [ $# -ne 3 ]; then
    echo "Usage: streamgc <genome.fasta> <coordinates.txt> <output.tsv>"
    exit 1
fi

genome_file=$1
coord_file=$2
output_file=$3

# Ensure samtools is installed
if ! command -v samtools &> /dev/null; then
    echo "Error: samtools is not installed or not in PATH."
    exit 1
fi

# Generate FASTA index if missing
if [ ! -e "${genome_file}.fai" ]; then
    samtools faidx $genome_file
fi

# Load chromosome lengths into an associative array
declare -A chr_lengths
while read -r chrom len _; do
    chr_lengths[$chrom]=$len
done < "${genome_file}.fai"

# Function to calculate GC content
calculate_gc() {
    local seq=$1
    local gc_count=$(echo "$seq" | tr -cd 'GCgc' | wc -c)
    local total_count=${#seq}
    if [ $total_count -eq 0 ]; then
        echo "0.0000"
    else
        echo "scale=2; $gc_count / $total_count" | bc
    fi
}

# Create header for output file
echo -e "Chromosome_Coordinate\tDistance\tGC_Ratio\tStrand" > $output_file

# Process each entry in the coordinate file
while IFS=$'\t' read -r chrom coord strand; do
    
    # Skip unknown chromosomes (not in genome index)
    chrom_len=${chr_lengths[$chrom]}
    if [ -z "$chrom_len" ]; then
        echo "Warning: Chromosome '$chrom' not found in genome index. Skipping..."
        continue
    fi

    # Loop through distances from -10,000 to +10,000 in steps of 100 bp
    for distance in $(seq -10000 100 10000); do
        
        # Calculate start and end positions of the window
        start=$((coord + distance))
        end=$((start + 99))
        
        # Skip windows that fall outside chromosome boundaries
        if (( start < 1 || end > chrom_len )); then
            continue
        fi
        
        # Extract sequence using samtools faidx and remove headers/newlines
        seq=$(samtools faidx "$genome_file" "${chrom}:${start}-${end}" | tail -n +2 | tr -d '\n')
        
        # Calculate GC content for the sequence and output results if valid sequence exists
        gc_ratio=$(calculate_gc "$seq")
        echo -e "${chrom}_${coord}\t${distance}\t${gc_ratio}\t${strand}" >> "$output_file"
        
    done
    
done < "$coord_file"

echo "GC content analysis complete. Results saved in $output_file"
