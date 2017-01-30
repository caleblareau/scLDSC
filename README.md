# scLDSC
LD Score Regression in Single Cells

You have to be careful with naming conventions else the the software will bark at you

```
python munge_sumstats.py --sumstats GIANT_BMI_Speliotes2010_publicrelease_HapMapCeuFreq.txt --merge-alleles w_hm3.snplist --out sumstats/BMI --a1-inc

python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.22 --ld-wind-cm 1 --print-snps list.txt --annot hemeBulk/heme22.annot --out hemeBulk/heme22

python ldsc.py --h2 sumstats/BMI.sumstats.gz --ref-ld-chr hemeBulk/heme --overlap-annot --w-ld-chr weights_hm3_no_hla/weights. --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --out BMI_baseline
```
