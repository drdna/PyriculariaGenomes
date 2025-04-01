# PyriculariaGenomes
Methods, Code and Data for the Project Genome Assemblies for Pyricularia species
## Secerted Protein Prediction
1. Run SignalP5 on maker protein models:
```bash
signalp -fasta <maker-proteins.fasta> -format short -prefix <genome-ID>
```
2. Count secreted proteins:
```bash
f=$(ls <signalP-summary-file>); echo ${f/_*/} | tr "\n" "\t"; awk '$2 ~ /^SP/' $f |  wc -l
```
## Variant Calling
SNPs were called from masked genome alignments using the StrictUnique4 module from iSNP caller.
1. BLAST masked reference genome against assemblies:
```bash
mkdir B71v2_BLAST
cd MASKED_GENOMES
for f in $(ls *fasta); do blastn -query ../B71v2_masked.fasta -subject $f -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > ../B71v2_BLAST/B71v2.${f/_*/}.BLAST; done
```
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
perl Generate_FASTA.pl PoTreeStrains.txt B71v2_SNPs B71v2.B71v2_BLAST_alignments
```
5. Build phylogenetic tree using RAXML:
```bash
code
```
