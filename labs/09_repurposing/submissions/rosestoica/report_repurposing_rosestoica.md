ex01
Fisier gasit: data/work/rosestoica/lab09/drug_gene_rosestoica.csv
Incarcat 71 interactiuni drug-gene
Numar medicamente: 40
    Abemaciclib: {'CDK4', 'CDK6'}
    Acalabrutinib: {'BTK'}
    Anastrozole: {'CYP19A1'}
Graf bipartit: 40 drugs, 32 genes, 71 edges
Graf salvat in: labs/09_repurposing/submissions/rosestoica/network_drug_gene_rosestoica.gpickle
Sumar medicamente salvat in: labs/09_repurposing/submissions/rosestoica/drug_summary_rosestoica.csv
      drug  num_targets
   Imatinib            4
   Olaparib            4
  Dasatinib            4
  Cisplatin            3
  Etoposide            3
  Rucaparib            3
  Sunitinib            3
  Sorafenib            3
Talazoparib            3
Doxorubicin            3

Muchii de similaritate (Jaccard >= 0.1): 44
Similaritati salvate in: labs/09_repurposing/submissions/rosestoica/drug_similarity_rosestoica.csv

Top 10 perechi de medicamente similare:
        drug1       drug2  similarity
  Abemaciclib Palbociclib         1.0
  Abemaciclib  Ribociclib         1.0
Acalabrutinib   Ibrutinib         1.0
  Anastrozole   Letrozole         1.0
   Bortezomib Carfilzomib         1.0
    Cetuximab   Erlotinib         1.0
    Cetuximab   Gefitinib         1.0
root@codespaces-4063eb:/workspaces/bioinf-y4-lab#  cd /workspaces/bioinf-y4-lab && python labs/09_repurposing/submissions/rosestoica/ex02_disease_proximity.py
 Fisier disease genes gasit: data/work/rosestoica/lab09/disease_genes_rosestoica.txt
 Graf incarcat din: labs/09_repurposing/submissions/rosestoica/network_drug_gene_rosestoica.gpickle
 Graf cu 40 medicamente si 32 gene
 Disease genes: 11 gene
    Gene: {'CDK6', 'CDK4', 'BRAF', 'BCL2', 'MDM2', 'BRCA2', 'ERBB2', 'BRCA1', 'TP53', 'EGFR', 'PARP1'}
 Disease genes prezente in graf: 11

 Calculam proximitatea medicamentelor fata de genele bolii...
 Ranking salvat in: labs/09_repurposing/submissions/rosestoica/drug_priority_rosestoica.csv

 Top 15 medicamente (proximitate maxima fata de disease genes):
       drug  distance  proximity_score
  Cisplatin  4.272727         0.228690
  Rucaparib  4.272727         0.228690
   Olaparib  4.272727         0.228690
Talazoparib  4.272727         0.228690
   Nutlin-3  4.636364         0.211132
  Etoposide  4.818182         0.203327
Doxorubicin  4.818182         0.203327
  Lapatinib  5.090909         0.192644
 Ribociclib  5.090909         0.192644
Palbociclib  5.090909         0.192644
Abemaciclib  5.090909         0.192644
 Pertuzumab  5.272727         0.186125
  Gefitinib  5.272727         0.186125
  Erlotinib  5.272727         0.186125
  Cetuximab  5.272727         0.186125

Bottom 5 medicamente (cele mai indepartate):
        drug  distance  proximity_score
Obinutuzumab       6.0         0.163934
Lenalidomide       6.0         0.163934
   Ibrutinib       6.0         0.163934
   Tamoxifen       6.0         0.163934
  Trametinib       6.0         0.163934


ex02
Fisier disease genes gasit: data/work/rosestoica/lab09/disease_genes_rosestoica.txt
Graf incarcat din: labs/09_repurposing/submissions/rosestoica/network_drug_gene_rosestoica.gpickle
Graf cu 40 medicamente si 32 gene
Disease genes: 11 gene
    Gene: {'CDK6', 'CDK4', 'BRAF', 'BCL2', 'MDM2', 'BRCA2', 'ERBB2', 'BRCA1', 'TP53', 'EGFR', 'PARP1'}
Disease genes prezente in graf: 11
  
Interpretare:

Medicamentele cu distanta mica (proximitate mare) au gene tinta care se suprapun cu genele bolii
PARP inhibitorii (Olaparib, Rucaparib, Talazoparib) sunt buni candidati deoarece tintesc BRCA1/2 si PARP1 - gene cheie in cancer
Medicamentele cu distanta mare (ex. Rituximab, Obinutuzumab) nu au suprapunere directa cu genele bolii