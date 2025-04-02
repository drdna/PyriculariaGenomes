# PyriculariaGenomes
Methods, Code and Data for the Project Genome Assemblies for Pyricularia species
## Trimming Sequence Reads
1. Quality was first assessed using FASTQC
2. Poor quality regions and adaptor contamination were then removed using Trimmomatic:
```bash
trimmomatic PE -phred33 -trimlog MyGenome_logfile.txt <MyGenome>_1.fq.gz <MyGenome>_2.fq.gz <MyGenome>_1_paired.fq <MyGenome>_1_unpaired.fq <MyGenome>_2_paired.fq <MyGenome>_2_unpaired.fq ILLUMINACLIP<path/to/adaptors.fasta>:2:30:10 SLIDINGWINDOW:20:20 MINLEN:120
```
## Genome Assembly
1. Genomes were assembled using a [velvetoptimiser](/scripts/velvetoptimiser_noclean.sh) script with k-mer values bracketing the value suggested by Velvet Advisor (https://dna.med.monash.edu/~torsten/velvet_advisor/) starting from k -40 to k +40. Initailly a step size of 10 was used:
```bash
sbatch velvetoptimser_noclean.sh <MyGenome_prefix> <starting_k> <ending_k> 10
```
2. After identifying the k value that produced the optimal assembly, this value was then bracketed from k - 10 to k +10, with a step size of 2:
```bash
sbatch velvetoptimser_noclean.sh <MyGenome_prefix> <starting_k> <ending_k> 2
```
3. Sequence headers were standardized usign the [SimpleFastaHeaders.pl](/scripts/SimpleFastaHeaders.pl) script:
```bash
perl SimpleFastaHeaders.pl contigs.fa <MyGenomeID>
```
4. And contigs < 200 nt in length were removed using [CullShortContigs.pl](/scripts/CullShortContigs.pl):
```bash
perl CullShortSequences.pl <MyGenomeID>_nh.fasta
```
## Genome Validation
1. Genome quality was assessed using BUSCO:
```bash
busco --in <MyGenome>_Final.fasta --out <MyGenome>_busco --mode genome --lineage_dataset ascomycota_odb10 -f
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
Generate metadata files for the _P. oryzae_ tree ([PoryzaeMetadata.txt](/data/PoryzaeMetadata.txt)) and the _Pyricularia_ tree ([PyriculariaMetadata.txt](/data/PyriculariaMetadata.txt)). Run the [Meyer_et_al_Fig1.R](/scripts/Meyer_et_al_Fig1.R) script on the support tree obtained from step 6.

