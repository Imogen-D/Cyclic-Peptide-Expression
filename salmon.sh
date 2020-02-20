#code to run Salmon on long read genomes with transcript and RNA seq data

salmon index -t <transcripts.fa> -i <strain>transcripts_index -k 31

#use k-mer as recommended. No decoy made
#single end reads

salmon quant -i <strain>transcripts_index -l A -r <strain>reads.fq --validateMappings -p 8 -o <strain>transcripts_quant

scp <strain>transcripts_quant/quant.sf ./<strain>quant.sf
#rename quantification file to be used in R script
