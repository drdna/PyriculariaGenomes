# PyriculariaGenomes

**Methods, Code and Data for the project Genome Assemblies for _Pyricularia_ species.
Steps 1 through 7 were completed by undergraduate students as independent research projects performed in the classes: ABT 480 / CS 485G: Special Topics - Applied Bioinformatics**
-------------------------------------------------------
## Table of Contents
1. [Assess sequence quality](#assess-sequence-quality)
2. [Trim poor quality sequence and adaptor contamination](#trim-poor-quality-sequence-and-adaptor-contamination)
3. [Genome assembly](#genome-assembly)
4. [Post processing of genome assembly for submission to NCBI](#post-processing-of-genome-assembly-for-submission-to-ncbi)
5. [Genome validation](#genome-validation)
6. [Identify genetic variants](#identify-genetic-variants)
7. [Gene prediction using MAKER](#gene-prediction-using-maker)
8. [Secreted protein prediction](#secreted-protein-prediction)
9. [Call haplotypes for each individual strain](#call-haplotypes-for-each-individual-strain)
10. [Generate a fasta file based on the haplotype calls](#generate-a-fasta-file-based-on-the-haplotype-calls)
11. [Generate phylogenetic trees](#generate-phylogenetic-trees)

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
1. Use Trimmomatic (v. 0.38) to remove poor quality regions and adaptor contamination. Note the following command includes parameters that perform well on most fungal datasets. The MINLEN parameter will need to be adjusted based on read lengths as well as fragment lengths as indicated by the distribution of adaptor contamination:
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
6. Run the [VelvetOptimiser](/scripts/velvetoptimiser_noclean.sh) script running VelvetOptimiser (v. 2.2.6) to perform assemblies at k-mer values bracketing the value suggested by  starting suggested by Velvet Advisor by 40 on each side (e.g. k -40 to k +40), with a step size of 10
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
3. Identify mitochondrial contigs by aligning to a reference mitochondrial genome. Save results for those contigs that align to mitochondrial sequences over > 90% of their length:
```
blastn -query MoMitochondrion.fasta -subject MyGenomeID_final.fasta -outfmt '6qseqid sseqid slen pident length mismatch gapopen qstart qend sstart send evalue score' | awk '$5/$3 > 0.9' > MoMitochondrion.MyGenomeID.BLAST
```
4. Create a list of mitochondrial contigs for uploading to the sequence assignment tab of the NCBI genome submission:
```bash
awk '{print $2 "\t" mitochondrion}' MoMitochondrion.MyGenomeID.BLAST > MyGenomeID_mitochondrial.csv
```
## Genome validation
1. Use the [BuscoSingularity.sh](/scripts/BuscoSingularity.sh) script to run BUSCO. The command line used is as follows: busco --in <MyGenome>_final.fasta --out <MyGenome>_busco --mode genome --lineage_dataset ascomycota_odb10 -f:
```bash
sbatch BuscoSingularity.sh <MyGenomeID>_final.fasta
```
## Gene prediction using MAKER
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
## Identify genetic variants
1. Use BLAST v. 2.16.0 to align \<MyGenomeID\>_final.fasta to a repeat-masked version of the B71 reference genome:
```
blastn -query B71v2sh_masked.fasta -subject <MyGenomeID>_final.fasta -evalue 1e-20 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' > B71v2sh.MyGenomeID.BLAST
```
2. Create a new directory for the BLAST results and copy the BLAST results into it:
```
mkdir MyGenomeID_BLAST
mv B71v2sh.MyGenomeID.BLAST MyGenomeID_BLAST
```
3. Copy the BLAST output file into the CLASS_BLASTS directory
4. Use the [CallVariants.sh](/scripts/CallVariants.sh) script to call the [StrictUnique.pm](/scripts/StrictUnique.pm) module from iSNPcaller to perform variant calling:
```
sbatch CallVariants.sh MyGenomeID_BLAST
```
5. Copy the SNP call output file into the CLASS SNPs directory

## Call haplotypes for each individual strain
1. Use the [Generate_haplotypes.pl](/scripts/Generate_haplotypes.pl) script to gather SNP calls for each genome and then check to make sure that SNPs is in a region that is unique in both the query genome and the reference. Print out two sets of haplotype data (1 = ref allele; 0 = alt allele; 9 = missing data/repeat region): a) all variant sites; and b) only sites called in every strain.
2. Arguments are as follows:
Generate_haplotypes.pl <strain-list> CLASS_SNPs CLASS_BLASTs <output_filename> <reference-sequence> <sequence-prefix>
```
perl Generate_haplotypes.pl PoABT480strains.txt CLASS_SNPS CLASS_BLASTS PoABT480Haplotypes B71v2sh_masked.fasta sequence
perl Generate_haplotypes.pl PyricABT480strains.txt CLASS_SNPS CLASS_BLASTS PyricABT480Haplotypes B71v2sh_masked.fasta sequence
```
## Generate a fasta file based on the haplotype calls
1. Use the Haplotypes2Fasta.pl script to generate the haplotype calls for a pre-determined list of strains:
```bash
perl Haplotypes2Fasta.pl PoTreeStrains.txt ClassHaplotypes.complete.txt
```
## Generate phylogenetic trees
1. Build a phylogenetic tree using RAxML v. 8.2.12:
```bash
fasta=<fasta_file>
prefix=<outputFilePrefix>
raxmlHPC-PTHREADS -T 12 -m GTRCAT -n $prefix -s $fasta -p 1234 -f a -x 4321 -# autoMRE
```
2. Add support values to nodes/branches using RAxML-NG v. 1.2.0:
```bash
besttree=<bestTreeFile>
bstrees=<bootstrapTreesFile>
prefix=<outputFilePrefix>
raxml-ng --support --tree $besttree --bs-trees $bstrees --threads 1 --prefix $prefix
```
7. Plot trees using ggtree and ggtreeextra:
Generate metadata files for the _P. oryzae_ tree ([PoryzaeMetadata.txt](/data/TreeBuilding/PoryzaeMetadata.txt)) and the _Pyricularia_ tree ([PyriculariaMetadata.txt](/data/TreeBuilding/PyriculariaMetadata.txt)). Open the files in FigTree and rotate branches to try and place all isolates from the same clade on branches that appear adjacent on the final tree. Run the [Meyer_et_al_Fig1.R](/scripts/Meyer_et_al_Fig1.R) script on the support trees obtained from step 6.

