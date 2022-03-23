###########
# Run duphold on output of gridss_post_process.py
# removes dels/dups based on fold change, 
# threshold reccomended by creators of duphold  
###########

# run duphold via slurm
VCF_PATH=target_vcf
BAM_PATH=target_bamfile
FASTA=reference_genome.fasta
OUTDIR=output_directory

VCF="$VCF_PATH""$SAMPLE"/"$SAMPLE"_simple_sv_gridss.vcf
BAM="$BAM_PATH""$SAMPLE".dedup.bam	
#run duphold
duphold -v "$VCF" -b "$BAM" -f "$FASTA" -t 1 -o "$SAMPLE"_duphold.vcf
#apply duphold filters
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3)' \
-o "$SAMPLE"_duphold_dupdel_filtered.vcf "$SAMPLE"_duphold.vcf
#get non dup/del
bcftools view -i '(SVTYPE = "INS") | (SVTYPE = "INV") | (SVTYPE = "BND")' -o "$SAMPLE"_ins_inv_bnd.vcf "$SAMPLE"_duphold.vcf

#get stats
DEL="$(grep "SVTYPE=DEL" "$VCF"|wc -l)"
DEL_filtered="$(grep "SVTYPE=DEL" "$SAMPLE"_duphold_dupdel_filtered.vcf|wc -l)" 
DUP="$(grep "SVTYPE=DUP" "$VCF"|wc -l)"
DUP_filtered="$(grep "SVTYPE=DUP" "$SAMPLE"_duphold_dupdel_filtered.vcf|wc -l)"
INV="$(grep "SVTYPE=INV" "$SAMPLE"_ins_inv_bnd.vcf|wc -l)"
INS="$(grep "SVTYPE=INS" "$SAMPLE"_ins_inv_bnd.vcf|wc -l)"
BND="$(grep "SVTYPE=BND" "$SAMPLE"_ins_inv_bnd.vcf|wc -l)"

# make filtered VCF - all the types
bgzip -c "$SAMPLE"_ins_inv_bnd.vcf>"$SAMPLE"_ins_inv_bnd.vcf.gz
bgzip -c "$SAMPLE"_duphold_dupdel_filtered.vcf>"$SAMPLE"_duphold_dupdel_filtered.vcf.gz
bcftools index "$SAMPLE"_ins_inv_bnd.vcf.gz
bcftools index "$SAMPLE"_duphold_dupdel_filtered.vcf.gz
bcftools concat -a "$SAMPLE"_ins_inv_bnd.vcf.gz "$SAMPLE"_duphold_dupdel_filtered.vcf.gz|bcftools sort -o "$SAMPLE"_filtered.vcf

#write out stats
TOTAL="$(grep -v "^#" "$SAMPLE"_filtered.vcf|wc -l)"
TOTAL_NOFILTER="$(grep -v "^#" "$VCF"|wc -l)"
echo sample total_no_fitler DEL DEL_filtered DUP DUP_filtered INV INS BND total> filter_stats.txt
echo "$SAMPLE" "$TOTAL_NOFILTER" "$DEL" "$DEL_filtered" "$DUP" "$DUP_filtered" "$INV" "$INS" "$BND" "$TOTAL">>filter_stats.txt
#remove temp files
rm "$SAMPLE"_ins_inv_bnd.vcf
rm "$SAMPLE"_ins_inv_bnd.vcf.gz
rm "$SAMPLE"_ins_inv_bnd.vcf.gz.csi
rm "$SAMPLE"_duphold_dupdel_filtered.vcf.gz
rm "$SAMPLE"_duphold_dupdel_filtered.vcf.gz.csi
rm "$SAMPLE"_duphold_dupdel_filtered.vcf
