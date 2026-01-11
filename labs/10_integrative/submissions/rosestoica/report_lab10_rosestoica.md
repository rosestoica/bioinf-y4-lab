EX01
python labs/10_integrative/submissions/rosestoica/ex01_PCA_and_viz.py
 SNP matrix: (50, 100) (features x samples)
 Expression matrix: (80, 100) (features x samples)
 Common samples: 100

 Normalizare z-score...
 Joint matrix: (130, 100) (features x samples)

 Rulare PCA...
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_snp_rosestoica.png
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_expression_rosestoica.png
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_joint_rosestoica.png

============================================================
COMPARATIE PCA - VARIANCE EXPLAINED
============================================================
Dataset                     PC1        PC2      Total
------------------------------------------------------------
SNP                       13.1%       5.1%      18.2%
Expression                11.0%       3.9%      14.9%
Joint (Integrated)        11.1%       3.3%      14.3%
===========================================================
Interpretare:

SNP-urile capteaza mai multa varianta (18.2%) in primele 2 componente
Integrarea dilueaza putin semnalul (14.3%) dar combina informatia din ambele straturi
PCA-ul integrat poate evidentia pattern-uri care nu sunt vizibile in date individuale

ex02
python labs/10_integrative/submissions/rosestoica/ex02_cross_omics.py
 SNP matrix: (50, 100)
 Expression matrix: (80, 100)
 Common samples: 100
 Calculare corelatii pentru 4000 perechi SNP-gene...
 Total perechi calculate: 4000

 Statistici corelatii:
    Mean |r|: 0.096
    Max |r|:  0.662
    Min |r|:  0.000

 Perechi cu |r| >= 0.3: 162
[OK] Salvat: labs/10_integrative/submissions/rosestoica/snp_gene_pairs_rosestoica.csv

 Top 10 perechi SNP-Gene corelate:
     snp     gene  correlation      p_value
rs100004 GENE_005     0.661676 6.656146e-14
rs100001 GENE_003     0.660357 7.761679e-14
rs100006 GENE_003     0.632127 1.746891e-12
rs100006 GENE_009     0.628437 2.564023e-12
rs100006 GENE_014     0.619697 6.236614e-12
rs100002 GENE_003     0.613088 1.199448e-11
rs100005 GENE_003     0.598281 4.920662e-11
rs100004 GENE_015     0.597152 5.463969e-11
rs100000 GENE_002     0.594029 7.283428e-11
rs100004 GENE_004     0.586824 1.397403e-10
[OK] Heatmap salvat: labs/10_integrative/submissions/rosestoica/correlation_heatmap_rosestoica.png

============================================================
SUMAR CROSS-OMICS ANALYSIS
============================================================
SNPs analizate:           50
Gene analizate:           80
Total perechi testate:    4000
Perechi semnificative:    162 (|r| >= 0.3)
SNPs cu corelatii:        14
Gene cu corelatii:        24
============================================================
root@codespaces-4063eb:/workspaces/bioinf-y4-lab#  cd /workspaces/bioinf-y4-lab && python labs/10_integrative/submissions/rosestoica/ex01_PCA_and_viz.py
 SNP matrix: (50, 100) (features x samples)
 Expression matrix: (80, 100) (features x samples)
 Common samples: 100

 Normalizare z-score...
 Joint matrix: (130, 100) (features x samples)
[OK] Salvat: labs/10_integrative/submissions/rosestoica/multiomics_concat_rosestoica.csv

 Rulare PCA...
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_snp_rosestoica.png
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_expression_rosestoica.png
[OK] Salvat: labs/10_integrative/submissions/rosestoica/pca_joint_rosestoica.png
COMPARATIE PCA - VARIANCE
============================================================
Dataset                     PC1        PC2      Total
------------------------------------------------------------
SNP                       13.1%       5.1%      18.2%
Expression                11.0%       3.9%      14.9%
Joint (Integrated)        11.1%       3.3%      14.3%


50 SNP-uri × 80 gene = 4000 perechi testate
162 perechi semnificative cu |r| ≥ 0.3
Corelatie maxima: r = 0.66 (rs100004 ↔ GENE_005)

Aceste corelatii arata ca anumite variante genetice (SNP-uri) sunt asociate cu niveluri modificate de expresie genica - exact ce cautam in studiile eQTL (expression Quantitative Trait Loci)!