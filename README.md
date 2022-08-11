# Kiwi-genomes

Commands, pipeline, and scripts for Genomic insights into the evolutionary relationships and demographic history of kiwi.

Note: This is not meant to be standalone commands but examples to be used for recreating the analyses in the manuscript

## Sequencing read filtering and mapping
- Script can be found here: https://github.com/Mvwestbury/Delphinoidea/blob/main/Modern_mapping_mem_PE.sh

## Finding sex-linked scaffolds/autosomes https://github.com/bioinfologics/satsuma2
- Use satsuma synteny

`satsuma - SatsumaSynteny -m 1 -n 30 -q Referencegenome.fasta -t Z_W_chromosomes.fasta -o output_directory`
- Obtain scaffold names

`cut -f 4 satsuma_summary.chained.out | sort | uniq | sed 's/1_/1\t/g' | cut -f 1 > Reference_XY.txt`
- Create text file with only autosomes and scaffolds >100kb

`samtools faidx Referencegenome.fasta`

`grep -v -f Reference_XY.txt Referencegenome.fasta.fai | awk '$2 >100000 {print $1}' > Reference_noXY_100kb.txt`

## Fasta call https://github.com/ANGSD/angsd

`angsd -dofasta 2 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -rf Reference_noXY_100kb.txt -i bamfile.bam -out bamfile_autosomes_100kb`

## Sliding window trees

- Generate fasta files for IQ-tree

-Create bed file with only scaffolds > 100kb
`awk '$2>99999 {print $1"\t1\t"$2}' Rowi.fa.fai > Genome.bed`

-Make bef file with 2kb windows with 1Mb slide
`bedtools makewindows -b Genome.bed -w 2000 -s 1000000 > Genome_w20kb_s1Mb.bed`

-Extract windows
`bedtools getfasta -fi Rowi.fa -bed Genome_w2kb_s1Mb.bed -fo Rowi_windows.fasta`

-Make directory for species of interest
`mkdir Rowi`

-Make each window into a unique fasta file
`cat Rowi_windows.fasta | awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > "Rowi/"filename}'`

-Add species name at the start of the fasta header
`sed -i 's/>/>Rowi_/g' Rowi/*`

- Repeat for all species

-Create a list with file names
`awk '$3>19999 {print $1":"$2"-"$3".fa"}' Genome_w20kb_s1Mb.bed > Filenames.txt`

-Combine single windows from all species into a single fasta file
`while read -r line; do cat Haastii/$line Mantelli/$line Owenii/$line Rowi/$line > Combined/Combined"_"$line; done < Filenames.txt`

- QuIBL

`QuIBL.py myInputFile2.txt`


## Dsuite https://github.com/millanek/Dsuite
`bcftools mpileup --threads 5 -f Referencegenome.fasta -Ou -R Scaffolds_100kb.txt -b bamlist.txt | bcftools call -mv -Ov --threads 5 -o All.vcf`
- Calculate the D (ABBA-BABA) and f4-ratio statistics for all possible trios of populations/species

`Dsuite Dtrios -t Tree.txt -o testDtri All_test.vcf Sets.txt`
- Calculate Fbranch

`Dsuite Fbranch Tree.txt testDtri_tree.txt > fbranch.txt`
- Plot Fbranch

`dtools.py fbranch.txt Tree.txt`

## Mutation rate https://github.com/ANGSD/angsd

`angsd -doIBS 2 -nthreads 5 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -minInd 9 -makematrix 1 -b bamlist.txt -rf Autosomes_100kb.txt -ref Referencegenome.fasta -out Kiwi_IBS` 

## hPSMC https://github.com/jacahill/hPSMC
- Merge 2x fasta consensuses into psmcfa file

`psmcfa_from_2_fastas.py -b10 -m5 species1_autosomes_100kb.fa species2_autosomes_100kb.fa > species1_species2_hPSMC.psmcfa`
- Run psmc

`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o species1_species2_hPSMC.psmc species1_species2_hPSMC.psmcfa`
- Plot PSMC to get predivergence Ne

`psmc_plot.pl -g 20 -u 7.90E-09 -R -s 10 test species1_species2_hPSMC.psmc`
- Run simulations

`hPSMC_quantify_split_time.py -N 20000 -l 500000 -u 5000000 -p 11 -s 11 -o output_prefix`


## Heterozygosity and runs of homozygosity https://github.com/grenaud/ROHan

`angsd -GL 1 -mininddepth 5 -minmapq 30 -minq 30 -uniqueonly 1 -only_proper_pairs 1 -docounts 1 -i bamfile -ref GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -P 5 -out species -doSaf 1 -anc GCF_000442215.1_Lipotes_vexillifer_v1_genomic.fna -rf Lipotes_vexillifer_noXY_100kb.txt -baq 1`


## PSMC https://github.com/lh3/psmc
- Create diploid consensus file

`samtools mpileup -Q 30 -uf reference.fasta bamfile.bam | bcftools call -c - | vcfutils.pl vcf2fq -d 10 | gzip > diploid.fq.gz`
- Create PSMC fasta file

`utils/fq2psmcfa -q20 diploid.fq.gz > diploid.psmcfa`
- Split PSMC fasta file for bootstrapping

`utils/splitfa diploid.psmcfa > diploid.split.psmcfa`
- Create PSMC data with whole sequence data

`psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o diploid.psmc diploid.psmcfa`
- Bootstrap replicates (with multiple threads)

`seq 100 | xargs -P 10 -i psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o diploid.round-{}.psmc diploid.split.psmcfa`
- Combine the psmc output files

`cat diploid.psmc diploid.round-*.psmc > combined.diploid.psmc`
