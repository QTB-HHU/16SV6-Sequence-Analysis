# 16SV6-Sequence-Analysis
Pipeline to process 16S (V6) sequence data coming from Ion-Torrent sequencer
Welcome to the 16SV6-Sequence-Analysis wiki!

The following steps describe the software and commands that were used in order to analyze the bacterial 16S V6 region obtain from a non-axenic culture of Phaeodactylum tricornutum.

Corresponding article: Dynamics of the bacterial community associated with Phaeodactylum tricornutum cultures Fiona Wanjiku Moejes, Ovidiu Popa, Antonella Succurro, Julie Maguireand Oliver Ebenhöh

Input:

Ion-Torrent FASTQ files containing raw reads without barcodes and adapter sequence. Primer are still attached and must be removed. Fasta file containing the primer pattern 'mapping' file containing sample information, needed for data evaluation.

Pipeline is adapted and modified from http://www.brmicrobiome.org/#!16s-profiling-ion-torrent-new/cpdg Data analysis for 16S microbial profiling from different benchtop sequencing platforms. J Microbiol Methods. 2014. doi: 10.1016/j.mimet.2014.08.018.

Software:

fastq-mcf Version: 1.04.807 (ref 1)

QIIMEfastaFormatter.pl (George Watts, University of Arizona 2013)

usearch v8.0.1517 (32Bit – opensource) (ref 2a,b)

Qiime v. 1.9.0 (ref 3)

Processing Steps:

Primer trimming
fastq-mcf -0 -f -l 20 -o $outfile primer.fa $input

Length trimming (based on the distribution of V6 from VAMPS (ref 4))
fastq-mcf -0 -l 52 -L 61 -o $outfile /dev/null $input

Quality filtering -USING USEARCH v8-
usearch -fastq_filter $PWD/reads2.fastq -fastq_maxee 0.5 -fastqout reads.fa

Dereplication -USING USEARCH v8-
usearch -derep_fulllength $PWD/reads.fa -output derep.fa -sizeout

Abundance sort and discard singletons -USING USEARCH v8-
usearch -sortbysize $PWD/derep.fa -output sorted.fa -minsize 2

OTU clustering using UPARSE method -USING USEARCH v8-
usearch -cluster_otus $PWD/sorted.fa -otus otus1.fa

Chimera filtering using reference rdp.gold.fa database -USING USEARCH v8-
usearch -uchime_ref $PWD/otus1.fa -db $PWD/rdp_gold.fa -strand plus -nonchimeras otus.fa

Map reads back to OTU database -USEARCH v8 script-
usearch -usearch_global $PWD/reads.fa -db $PWD/otus.fa -strand plus -id 0.97 -uc map.uc

Assign taxonomy to OTUS using uclust method on QIIME (use the file “otus.fa” from UPARSE as input file)
to speed up use parallel_assign_taxonomy_uclust.py from QIIME

assign_taxonomy.py -i $PWD/otus.fa -o output -r $PWD/rep_set/97_otus.fasta -t $PWD/taxonomy/97_otu_taxonomy.txt

Align sequences on QIIME, using SILVA reference sequences (use the file “otus.fa” from UPARSE as input file)
align_seqs.py -i $PWD/otus.fa -o rep_set_align -t $PWD/rep_set_aligned/97_otus.fasta

Filter alignments on QIIME
filter_alignment.py -i $PWD/otus_aligned.fasta -o filtered_alignment

Make the reference tree on QIIME
make_phylogeny.py -i $PWD/otus_aligned_pfiltered.fasta -o rep_set.tre

Convert UC to otu-table.txt -UPARSE PYTHON SCRIPT-
python $PWD/uc2otutab.py $PWD/map.uc > otu_table.txt

Convert otu_table.txt to otu-table.biom, used by QIIME -BIOM SCRIPT-
biom convert -i $PWD/otu_table.txt -o otu_table.biom --table-type="otu table"

Add metadata (taxonomy) to OTU table
biom add-metadata -i $PWD/otu_table.biom -o otu_table_tax.biom --observation-metadata-fp $PWD/otus_tax_assignments.txt --observation-header OTUID,taxonomy,confidence --sc-separated taxonomy --float-fields confidence

Check OTU Table on QIIME
biom summarize-table -i $PWD/otu_table_tax.biom -o results_biom_table

Run diversity analyses on QIIME (or any other analysis of your choice). The parameter “-e” is the sequencing depth to use for even sub-sampling and maximum rarefaction depth. You should review the output of the ‘biom summarize-table’ (step 15) command to decide on this value.
core_diversity_analyses.py -i $PWD/otu_table_tax.biom -m $PWD/mapping_file.txt -t $PWD/rep_set.tre -e xxxx -o $PWD/core_output

References

Erik Aronesty (2011). ea-utils : "Command-line tools for processing biological sequencing data"; http://code.google.com/p/ea-utils 2a) Edgar,RC (2010) Search and clustering orders of magnitude faster than BLAST, Bioinformatics 26(19), 2460-2461. doi: 10.1093/bioinformatics/btq461 http://drive5.com/usearch 2b) Edgar, R.C. (2013) UPARSE: Highly accurate OTU sequences from microbial amplicon reads, Nature Methods [Pubmed:23955772, dx.doi.org/10.1038/nmeth.2604] http://drive5.com/usearch
J Gregory Caporaso, Justin Kuczynski, Jesse Stombaugh, Kyle Bittinger, Frederic D Bushman, Elizabeth K Costello, Noah Fierer, Antonio Gonzalez Pena, Julia K Goodrich, Jeffrey I Gordon, Gavin A Huttley, Scott T Kelley, Dan Knights, Jeremy E Koenig, Ruth E Ley, Catherine A Lozupone, Daniel McDonald, Brian D Muegge, Meg Pirrung, Jens Reeder, Joel R Sevinsky, Peter J Turnbaugh, William A Walters, Jeremy Widmann, Tanya Yatsunenko, Jesse Zaneveld and Rob Knight; Nature Methods, 2010; doi:10.1038/nmeth.f.303 http://qiime.org/index.html
Huse SM, Mark Welch DB, Voorhis A, Shipunova A, Morrison HG, et al. (2014) VAMPS: a website for visualization and analysis of microbial population structures. BMC Bioinformatics 15: 41. https://vamps.mbl.edu/resources/databases.php
