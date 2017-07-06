# scLDSC
LD Score Regression in Single Cells

You have to be careful with naming conventions else the the software will bark at you

```
python munge_sumstats.py --sumstats sumstats/cd.u.txt --N 27726 --merge-alleles sumstats/w_hm3.snplist --out sumstats/CROHN --a1-inc

python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.22 --ld-wind-cm 1 --print-snps list.txt --annot hemeBulk/heme22.annot --out hemeBulk/heme22

python ldsc.py --h2 sumstats/BMI.sumstats.gz --ref-ld-chr hemeBulk/heme --overlap-annot --w-ld-chr weights_hm3_no_hla/weights. --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out BMI_baseline
```


## Single Cells:
```
bsub -q big -M 64000 -R 'rusage[mem=64000]' python ldsc.py --h2 sumstats/BMI.sumstats.gz --ref-ld-chr scMyeloid/scHeme --overlap-annot --w-ld-chr weights_hm3_no_hla/weights. --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out scMyeloid/BMI
```

## What's going on 

The [raw.annot](raw.annot) folder contains annotation files for 1000 Genomes Project. Here's what
the first few lines look like--

```
zcat < ../raw.annot/raw.1.annot.gz | head -5
CHR	BP	SNP	CM	base
1	11008	rs575272151	0.0	1
1	11012	rs544419019	0.0	1
1	13110	rs540538026	0.0	1
1	13116	rs62635286	0.0	1
```

What one would want to do for another set of data is take the annotation files and use `awk`
to make similar files and then plug and chug the key Rscript `code/bed2annotate.R`

