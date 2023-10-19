
# Build singularity image
singularity build --fakeroot ParTIES.sif ParTIES.def

# create parties alias
alias parties='singularity exec  ParTIES.sif parties'

# test ParTIES
EXEMPLE_DIR=../example/
gunzip $EXEMPLE_DIR/*fastq.gz
	
FASTQ1=$EXEMPLE_DIR/Example_reads_1.fastq
FASTQ2=$EXEMPLE_DIR/Example_reads_2.fastq
OUT=Test2
	
# Mapping using BOWTIE2
parties Map -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT \
	 -fastq1 $FASTQ1 -fastq2 $FASTQ2 \
	 -max_insert_size 500 -index_genome -threads 4  -v

# MIRAA detection
parties MIRAA -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT \
	 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
	 -min_break_coverage 5 -threads 4  -v

parties Assembly -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT \
	 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
	 -miraa $OUT/MIRAA/MIRAA.gff3 \
	 -fastq1 $FASTQ1 -fastq2 $FASTQ2 \
	 -insert_size 300 -kmer 51 -threads 4  -v
 
parties MICA -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT  -skip_repeat_masker  \
	 -bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
	 -miraa $OUT/MIRAA/MIRAA.gff3 \
	 -germline_genome $OUT/Assembly/VELVET_51_at_least_one_no_match/VELVET_51_at_least_one_no_match_contigs.fa \
	 -germline_genome $OUT/Assembly/VELVET_51_no_filter/VELVET_51_no_filter_contigs.fa \
	 -germline_genome $OUT/Assembly/VELVET_51_no_mac_junctions/VELVET_51_no_mac_junctions_contigs.fa \
	 -insert_size 300 -threads 4  -v


parties Insert -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT \
	-ies $OUT/MICA/MICA.gff3 -suffix _with_IES \
	-threads 4 -v 

parties Map -genome $OUT/Insert/Insert.fa -out_dir $OUT \
	-fastq1 $FASTQ1 -fastq2 $FASTQ2 \
	-max_insert_size 500 -index_genome -threads 4 -v -force



parties MIRET -genome $EXEMPLE_DIR/scaffold51_1.fa -out_dir $OUT \
	-germline_genome $OUT/Insert/Insert.fa \
	-bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
	-germline_bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
	-ies $OUT/MICA/MICA.gff3 \
	-germline_ies $OUT/Insert/Insert.gff3 \
	-score_method Boundaries -threads 4 -v


parties MILORD -genome $OUT/Insert/Insert.fa -out_dir $OUT \
	-bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
	-ies $OUT/Insert/Insert.gff3 \
	-threads 4 -v 


parties Compare -genome $OUT/Insert/Insert.fa -out_dir $OUT \
	-reference_set $OUT/Insert/Insert.gff3 \
	-current_set $OUT/MILORD/MILORD.gff3 \
	-threads 4 -v 


parties Concatemer -genome $OUT/Insert/Insert.fa -out_dir $OUT \
	-seq_id null -fastq1 $FASTQ1 -fastq2 $FASTQ2 \
	-ies $OUT/MICA/MICA.gff3 \
	-bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam \
	-germline_bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
	-threads 4 -v


parties MEND -out_dir $OUT \
	-bam $OUT/Map/$OUT.scaffold51_1.fa.BOWTIE.sorted.bam  \
	-genome $EXEMPLE_DIR/scaffold51_1.fa \
	-germline_bam $OUT/Map/$OUT.Insert.fa.BOWTIE.sorted.bam \
	-germline_genome $OUT/Insert/Insert.fa \
	-germline_ies $OUT/Insert/Insert.gff3 \
	-ies $OUT/MICA/MICA.gff3 \
	-excision_errors $OUT/Compare/Compare.current.gff3 \
	-threads 4 -v 
 
 
