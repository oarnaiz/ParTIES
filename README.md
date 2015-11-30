# ParTIES
PARamecium Toolbox for Interspersed DNA Elimination Studies

# Description
We present the Paramecium Toolbox for Interspersed DNA Elimination Studies (ParTIES), designed for Paramecium species, 
that (i) identifies eliminated sequences, (ii) measures their presence in a sequencing sample and (iii) detects rare elimination polymorphisms.

For a full description of the software, its options and results look at the user manual ("user_manual.pdf")





# Install
ParTIES requires some other programs. An installation example is provided at the end of the user manual, as well as in the "INSTALL" file. 
Once all the dependencies are installed, use the "check" file to make sure everything is ok.
```bash
 ./check
```





# Usage
To run ParTIES use the following command
```bash
parties [MODE] : PARamecium Toolbox for Interspersed DNA Elimination Studies
 Run : Run ParTIES using the configuration file
 Map : Map reads on a reference using bowtie2
 MIRAA : Method of Identification by Read Alignment Anomalies
 MICA : Method of Identification by Comparison of Assemblies
 Insert : Insert IES within a genome to create an IES containing reference
 MIRET : Method of Ies RETention
 Assembly : Filter reads and assemble them
 MILORD : Method of Identification and Localization of Rare Deletion
 Compare : Compare IES/InDel datasets
```





# License
This software is distributed under the GNU GPL v3 license. See the "LICENSE" file for details.




# How to cite
If you use this software please cite the following publication

Denby Wilkes C, Arnaiz O, Sperling L. ParTIES : a toolbox for Paramecium interspersed DNA elimination studies. 
Bioinformatics. 2015 Nov 20. pii: btv691. [Epub ahead of print] PubMed PMID: 26589276.



# Example
The example directory contains the following files :

| File | Description |
| ----------- | ----------- |
| scaffold51_1.fa | The somatic reference for this example [a single somatic scaffold of the Paramecium tetraurelia genome] |
| Example_reads_1.fastq.gz | Read file 1 (paired with read file 2), 100 nt-long reads |
| Example_reads_2.fastq.gz | Read file 2 (paired with read file 1), 100 nt-long reads |
| example.cfg | Configuration file used to gives all needed options to run PARTIES with the All module (see below) |


The following command will run the entire pipeline (based on the config file), generating results in the $OUT  directory. You can comment lines in the configuration file by adding "#" at the begining of a line.


```bash
OUT=Test1
gunzip example/Example_reads_1.fastq.gz
gunzip example/Example_reads_2.fastq.gz
parties Run -genome example/scaffold51_1.fa -out_dir $OUT -config example/example.cfg
```

Before running the pipeline, check the configuration file ("example/example.cfg") to set the number of threads.

You can also run each step independently, specifying the intermediate result files on the command line.

The map module will align the reads on the reference.
```bash
OUT=Test2
parties Map -genome example/scaffold51_1.fa -out_dir $OUT \
 -fastq1 example/Example_reads_1.fastq -fastq2 example/Example_reads_2.fastq \
 -max_insert_size 500 -index_genome -threads 4 
```

The MIRAA module searches for breakpoints in an alignment file.
```bash
parties MIRAA -genome example/scaffold51_1.fa -out_dir $OUT \
 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
 -min_break_coverage 5 -threads 4 
```

The Assembly module filters sequencing reads and assemble them into contigs. Three different assemblies are created.
```bash
parties Assembly -genome example/scaffold51_1.fa -out_dir $OUT \
 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
 -miraa $OUT/MIRAA/MIRAA.gff3 \
 -fastq1 example/Example_reads_1.fastq -fastq2 example/Example_reads_2.fastq \
 -insert_size 300 -kmer 51 -threads 4 
```

The MICA module computes comparisons between genomes, looking for insertions in germline genomes.
```bash
parties MICA -genome example/scaffold51_1.fa -out_dir $OUT \
 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
 -miraa $OUT/MIRAA/MIRAA.gff3 \
 -germline_genome $OUT/Assembly/VELVET_51_at_least_one_no_match/VELVET_51_at_least_one_no_match_contigs.fa \
 -germline_genome $OUT/Assembly/VELVET_51_no_filter/VELVET_51_no_filter_contigs.fa \
 -germline_genome $OUT/Assembly/VELVET_51_no_mac_junctions/VELVET_51_no_mac_junctions_contigs.fa \
 -insert_size 300 -threads 4 
```


The Insert module creates an IES containing reference
```bash
parties Insert -genome example/scaffold51_1.fa -out_dir $OUT \
 -ies $OUT/MICA/MICA.gff3 -suffix _with_IES \
 -threads 4 
```

We use the Map module once again to align the reads on the IES containing reference.
```bash
parties Map -genome $OUT/Insert/Insert.fa -out_dir $OUT \
 -fastq1 example/Example_reads_1.fastq -fastq2 example/Example_reads_2.fastq \
 -max_insert_size 500 -index_genome -threads 4 -force
```

The MIRET module calculates precisely the level of retention of each IES in a sample.
```bash
parties MIRET -genome example/scaffold51_1.fa -out_dir $OUT \
 -germline_genome $OUT/Insert/Insert.fa \
 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
 -germline_bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
 -ies $OUT/MICA/MICA.gff3 \
 -germline_ies $OUT/Insert/Insert.gff3 \
 -score_method Boundaries -threads 4 
```

The MILORD module searches for rare deletions in sequencing reads compared to a reference.

When run on a germline genome, we do expect to see deletions that correspond to somatic reads.
```bash
parties MILORD -genome $OUT/Insert/Insert.fa -out_dir $OUT \
 -bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
 -ies $OUT/Insert/Insert.gff3 \
 -threads 4 
```

The Compare module allows coordinate-based comparisons between elements (MICA and/or MILORD results)
```bash
parties Compare -genome $OUT/Insert/Insert.fa -out_dir $OUT \
 -reference_set $OUT/Insert/Insert.gff3 \
 -current_set $OUT/MILORD/MILORD.gff3 \
 -threads 4 
```

