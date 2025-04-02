# PyriculariaGenomes
Methods, Code and Data for the Project Genome Assemblies for Pyricularia species
## Trimming Sequence Reads
Quality was first assessed using FASTQC and poor quality regions and adaptor contamination were then removed using Trimmomatic:
```bash
trimmomatic PE -phred33 -trimlog MyGenome_logfile.txt <MyGenome>_1.fq.gz <MyGenome>_2.fq.gz <MyGenome>_1_paired.fq <MyGenome>_1_unpaired.fq <MyGenome>_2_paired.fq <MyGenome>_2_unpaired.fq ILLUMINACLIP<path/to/adaptors.fasta>:2:30:10 SLIDINGWINDOW:20:20 MINLEN:120
```
## Genome Assembly
Genomes were assembled using a [velvetoptimiser](/scripts/velvetoptimiser_noclean.sh) script with k-mer values ranging around the value suggested by Velvet Advisor (https://dna.med.monash.edu/~torsten/velvet_advisor/):
```bash
sbatch velvetoptimser_noclean.sh <MyGenome_prefix> <starting_k> <ending_k> <stepsize>
```
## Secreted Protein Prediction
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
fasta=<fasta_file>
prefix=<outputFilePrefix>
raxmlHPC-PTHREADS /Applications/standard-RAxML-master/raxmlHPC-PTHREADS -T 12 -m GTRCAT -n $prefix -s $fasta -p 1234 -f a -x 4321 -# autoMRE
```
6. Add support values to nodes/branches:
```bash
besttree=<bestTreeFile>
bstrees=<bootstrapTreesFile>
prefix=<outputFilePrefix>
raxml-ng --support --tree $besttree --bs-trees $bstrees --threads 1 --prefix $prefix
```
7. Plot trees using ggtree and ggtreeextra:
Generate metadata files for the _P. oryzae_ tree ([PoryzaeMetadata.txt](/data/PoryzaeMetadata.txt)) and the _Pyricularia_ tree ([PyriculariaMetadata.txt](/data/PyriculariaMetadata.txt)). Run the [Meyer_et_al_Fig1.R](/scripts/Meyer_et_al_Fig1.R) script on the support tree obtained from step 6.

