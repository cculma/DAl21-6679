# DAl21-6679

Project with DArT data of GSDIG for GWAS

## Pipeline to calculate the LD decay

I am working with DAl21_AMCC.vcf (1502 gen by 2953 markers) instead of DAl21_MADC.vcf (1479 gen by 6494 SNPs). $M^{(n \times m)}$ where $n$ is the number of genotypes (accessions) and $m$ is the number of SNPs.
1. VCF to updog produce an multi-updog element: `vcf_updog_2.R`
2. multi-updog element is converted to a ldsep element to calculate the $R^2$ and then the ld decay: `LD_plot.R` and/or `LD_plot1.R`.
- ldsep element is a matrix $m \times m$ with $R^2$ values. We are only interested in compare the $R^2$ among chomosomes (linkage groups). This correspond to lines 72 to 120 of `LD_plot.R`.
