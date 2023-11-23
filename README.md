# PyriculariaGenomes
Methods, Code and Data for the Project Genome Assemblies for Pyricularia species
## Variant Calling
SNPs were called from masked genome alignments using the StrictUnique4 module from iSNP caller.
1. BLAST masked reference genome against assemblies:
```bash
mkdir B71v2_BLAST
cd MASKED_GENOMES
for f in `ls *fasta | grep -v ^Cr`; do blastn -query ../B71v2_masked.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../B71v2_BLAST/B71v2.${f/_*/}.BLAST
for f in `ls Cr*fasta`; do blastn -query ../B71v2_masked.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../B71v2_BLAST/B71v2.${f/_*/}.BLAST
2. Call SNPs:
```bash
cd ..
perl Run_SU4.pl B71v2_BLAST B71v2_SNPs
```
3. Mask sites in reference genome that were not queried in ALL genomes:
```bash
perl Create_alignment_strings_multiBLAST.pl B71v2_masked.fasta B71v2_BLAST
```
4. Generate .fasta file from SNP calls:
```bash
perl Generate_FASTA.pl Meyer_et_al_StrainList.txt B71v2_SNPs B71v2.B71v2_BLAST_alignments
```
5. Build phylogenetic tree using RAXML:
