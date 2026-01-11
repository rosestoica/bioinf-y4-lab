=== Gene Co-Expression Network Analysis ===
Dataset: GSE48350 de la NCBI GEO
Output: /workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/rosestoica/modules_rosestoica.csv

0. Descărcare date de la NCBI GEO...
   Descărcare de la: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE48nnn/GSE48350/matrix/GSE48350_series_matrix.txt.gz
   Salvat în: /workspaces/bioinf-y4-lab/data/work/rosestoica/lab06/expression_matrix.csv

1. Citire matrice de expresie...
   Dimensiuni inițiale: 54675 gene × 253 probe

2. Preprocesare (log2 + filtrare varianță + top genes)... prag=1.0
   Gene rămase după filtrare varianță: 109 din 54675

3. Calculare matrice de corelație (spearman)...
   Dimensiuni matrice corelație: (109, 109)

4. Construire matrice de adiacență (prag=0.7)...
   Număr muchii potențiale: 246

5. Construire graf NetworkX...
   Graf creat cu 75 noduri și 246 muchii.

6. Detectare module (comunități)...
Folosim algoritmul Louvain pentru detectarea modulelor.
   S-au detectat 11 module.

7. Salvare rezultate...
   Am salvat mapping-ul gene→modul în: /workspaces/bioinf-y4-lab/labs/06_wgcna/submissions/rosestoica/modules_rosestoica.csv

=== Analiză completă ===


Rețeaua păstrează topologia (cine este hub, distribuția gradelor, conexiuni indirecte), permite măsuri de rețea (centralitate, betweenness, conectivitate). - arata relatiile/interactiunile dintre gene

Clustering reduce totul la atribuirea fiecărei gene unui cluster, pierzând detalii despre relațiile locale dintre gene. -  arata segmenatarea pe categorii