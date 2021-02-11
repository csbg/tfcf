# Read matching sequence (RNA-seq)
gunzip -c SLX-19954.SITTH12.HY75WDRXX.s_1.r_2.fq.gz | grep "TATTTCTAGCTCTAAAAC" | sort > x_TATTTCTAGCTCTAAAAC.fq
gunzip -c SLX-19954.SITTH12.HY75WDRXX.s_1.r_2.fq.gz | grep "ATAAAGATCGAGATTTTG" | sort > x_ATAAAGATCGAGATTTTG.fq

# Read matching sequence (guide reads)
gunzip -c SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz | grep "CTCTAAAAC" | sort > x2_CTCTAAAAC.fq
gunzip -c SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz | grep "CTCTTAAAC" | sort > x2_CTCTTAAAC.fq
gunzip -c SLX-19954.HY75WDRXX.s_1.r_2.lostreads.fq.gz | grep "TATTTCTAGCTCTAAAAC" | sort > x2_TATTTCTAGCTCTAAAAC.fq

Rscript ANALYSIS.R

# Top 5000 reads sorted
gunzip -c SLX-19954.SITTH12.HY75WDRXX.s_1.r_2.fq.gz | head -20000 | perl -ne 'print unless (0 != ($.-2) % 4)' | sort > x2.fq

# FastQC
gunzip -c SLX-19954.SITTH12.HY75WDRXX.s_1.r_2.fq.gz | head -50000 > test.fq
gzip test.fq
/Applications/FastQC.app/Contents/MacOS/fastqc test.fq.gz