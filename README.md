# PyriculariaGenomes
Methods, Code and Data for the Project Genome Assemblies for Pyricularia species
## Table of Contents
1. [Assess sequence quality](#assess-sequence-quality)
2. [Trim poor quality sequence and adaptor contamination](#trim-poor-quality-sequence-and-adaptor-contamination)
3. [Post processign of genome assembly for submission to NCBI](#post-processing-of-genome-assembly-for-submission-to-ncbi)
4. [Identify genetic variants](#identify-genetic-variants)

## Assess Sequence Quality
1. Log onto VM
2. Change into *sequences* directory
3. Run FASTQC (v. 0.9.11) on sequence reads
```
fastqc MyGenomeID*fq.gz
```
4. Secure copy the resulting .html files to your local machine and double click to view reports in your default web browser
5. Inspect FASTQC reports to look for poor quality sequence and adaptor contamination at the ends of reads
6. Decide on suitable parameters for trimming (depends on individual sequence dataset)

## Trim poor quality sequence and adaptor contamination
1. Use Trimmomatic to remove poor quality regions and adaptor contamination. Note the following command includes parameters that perform well on most fungal datasets. The MINLEN parameter will need to be adjusted based on read lengths as well as fragment lengths as indicated by the distribution of adaptor contamination:
```bash
trimmomatic PE -phred33 -trimlog MyGenome_logfile.txt <MyGenome>_1.fq.gz <MyGenome>_2.fq.gz <MyGenome>_1_paired.fq <MyGenome>_1_unpaired.fq <MyGenome>_2_paired.fq <MyGenome>_2_unpaired.fq ILLUMINACLIP<path/to/adaptors.fasta>:2:30:10 SLIDINGWINDOW:20:20 MINLEN:120
```
## Genome Assembly
1. Transfer trimmed reads to the Morgan Compute Cluster using **scp**
2. Create a designated assembly directory
```
mkdir MyGenomeID
```
3. Copy the reads into the assembly directory
```
cp MyGenomeID*_paired.fq.gz MyGenomeID
```
4. Use Velvet Advisor (https://dna.med.monash.edu/~torsten/velvet_advisor/) to determine a reasonable k-mer value based on metrics for your trimmed reads
5. Change into the MyGenomeID directory
6. Run the [VelvetOptimiser](/scripts/velvetoptimiser_noclean.sh) script running VelvetOptimiser v. 2.2.6 to perform assemblies at k-mer values bracketing the value suggested by  starting suggested by Velvet Advisor by 40 on each side (e.g. k -40 to k +40), with a step size of 10
```bash
sbatch velvetoptimser_noclean.sh <MyGenome_prefix> <starting_k> <ending_k> 10
```
7. Inspect the log file (XX-XX-XXXX-YY-YY-YY_Logfile.txt, where the Xs and Ys correspond to values for date and time, respectively). Identify the k value that produced the optimal assembly and perform a new assembly bracketing this value from optimalK - 10 to optimalK +10, using a step size of 2:
```bash
sbatch velvetoptimser_noclean.sh <MyGenome_prefix> <starting_k> <ending_k> 2
```
## Post processing of genome assembly for submission to NCBI
1. Remove any contigs less than 200 nt in length using the [CullShortContigs.pl](/scripts/CullShortContigs.pl) script:
```bash
perl CullShortSequences.pl <MyGenomeID>.fasta
```
2. Re-number contigs for easier downstream processing using the [SimpleFastaHeaders.pl](/scripts/SimpleFastaHeaders.pl) script:
```
perl SimpleFastaHeaders.pl <MyGenomeID>_temp.fasta
```
3. Identify mitochondrial contigs by aligning to a reference mitochondrial genome:
```
blastn -query MoMitochondrion.fasta -subject MyGenomeID_final.fasta -outfmt '6qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue score' | awk '$5/$3 > 0.9' > MoMitochondrion.MyGenomeID.BLAST
```
## Identify genetic variants
1. Use BLAST v. 2.16.0 to align <MyGenome>_final.fasta to a repeat-masked version of the B71 reference genome:
```
blastn -query B71v2sh_masked.fasta -subject <MyGenomeID>_final.fasta -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > B71v2sh.MyGenomeID.BLAST
```
2. Create a new directory for the BLAST results and copy the BLAST results into it:
```
mkdir MyGenomeID_BLAST
mv B71v2sh.MyGenomeID.BLAST MyGenomeID_BLAST
```
3. Use the [CallVariants.sh](/scripts/CallVariants.sh) script to call the [StrictUnique.pm](StrictUnique.pm) module from iSNPcaller to perform variant calling:
```
sbatch CallVariants.sh MyGenomeID_BLAST
```
   
## Genome Validation
1. Use the [BuscoSingularity.sh](BuscoSingularity.sh) script to run BUSCO. The command line used is as follows: busco --in <MyGenome>_final.fasta --out <MyGenome>_busco --mode genome --lineage_dataset ascomycota_odb10 -f:
```bash
sbatch BuscoSingularity.sh <MyGenomeID>_final.fasta
```
## Gene Prediction using MAKER
1. The MAKER configuration files were created:
```bash
maker -CTL
```
2. The following settings were modified in the maker_opts.ctl file:
- **genome** = `/home/<username>/genes/<MyGenome>_final.fasta`
- **model_org** = *(set to blank)*
- **repeat_protein** = *(set to blank)*
- **snaphmm** = `/home/<username>/genes/Moryzae.hmm`
- **augustus_species** = `magnaporthe_grisea`
- **keep_preds** = `1`
- **protein** = `/home/yourusername/genes/genbank/ncbi-protein-Magnaporthe_organism.fasta`

3. Then run maker:
```bash
maker 2>&1 | tee maker.log
```
## Secreted Protein Prediction
1. Run SignalP5 on maker protein models:
```bash
fasta_merge -d <MyGenomeID>.maker.output/<MyGenomeID>_master_datastore_index.log -o <MyGenomeID>_genes
signalp -fasta <MyGenomeID>_genes/<MyGenomeID>-proteins.fasta> -format short -prefix <MyGenomeID>
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
Generate metadata files for the _P. oryzae_ tree ([PoryzaeMetadata.txt](/data/TreeBuilding/PoryzaeMetadata.txt)) and the _Pyricularia_ tree ([PyriculariaMetadata.txt](/data/TreeBuilding/PyriculariaMetadata.txt)). Open the files in FigTree and rotate branches to try and place all isolates from the same clade on branches that appear adjacent on the final tree. Run the [Meyer_et_al_Fig1.R](/scripts/Meyer_et_al_Fig1.R) script on the support trees obtained from step 6.

