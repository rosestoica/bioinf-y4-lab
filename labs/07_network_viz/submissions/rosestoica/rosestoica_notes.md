# Lab 07 - Network Visualization Notes

Metoda de layout
Am folosit spring layout (algoritmul Fruchterman-Reingold):

pos = nx.spring_layout(G, seed=SEED, k=2/np.sqrt(G.number_of_nodes()))

Ce avantaje aduce vizualizarea fata de analiza numerica din Lab 6?

1. **Perspectiva globala asupra retelei** - Reprezentarea grafica permite identificarea imediata a grupurilor de gene si a conexiunilor dintre ele, mult mai repede deact reprezantarea tabelara cu valori

2. **Evidentierea hub-urilor** - Genele cu grad mare de conectivitate ies in evidenta prin dimensiune si pozitie centrala, oferind o selectie naturala a candidatilor importanti.

3. **Descoperirea structurilor ascunse** - Aranjamentul spatial dezvaluie clustere si relatii intre module care raman invizibile in analiza tabelara.


---

## Output executie

```
 python labs/07_network_viz/submissions/ex_01_network_viz.py
=== Network Visualization & Hub Genes ===
Handle: rosestoica

1. Verificare fisiere input...
   ✓ Gasit: expression_matrix.csv
   ✓ Gasit: modules_rosestoica.csv

2. Incarcare date...
   Matrice expresie: 54675 gene × 253 probe
   Module incarcate: 75 gene in 11 module

3. Preprocesare si construire adiacenta...
   Gene dupa filtrare: 109
   Adiacenta calculata: (109, 109)

4. Construire graf...
   Graf: 75 noduri, 246 muchii

5. Asignare culori dupa modul...
/workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/ex_01_network_viz.py:134: MatplotlibDeprecationWarning: The get_cmap function was deprecated in Matplotlib 3.7 and will be removed in 3.11. Use ``matplotlib.colormaps[name]`` or ``matplotlib.colormaps.get_cmap()`` or ``pyplot.get_cmap()`` instead.
  cmap = plt.cm.get_cmap('tab10')
   Noduri colorate: 75/75

6. Calculare hub genes (top 10)...
   Hub genes:
      202376_at: degree=18, betweenness=0.0711
      204489_s_at: degree=18, betweenness=0.0457
      212063_at: degree=17, betweenness=0.0435
      209835_x_at: degree=17, betweenness=0.0388
      208651_x_at: degree=15, betweenness=0.0338
      238123_at: degree=14, betweenness=0.0473
      204490_s_at: degree=14, betweenness=0.0327
      229797_at: degree=14, betweenness=0.0211
      227253_at: degree=14, betweenness=0.0095
      1557905_s_at: degree=13, betweenness=0.0141

7. Vizualizare retea...

8. Salvare figura...
   Salvat: /workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/network_rosestoica.png

9. Salvare hub genes...
   Salvat: /workspaces/bioinf-y4-lab/labs/07_network_viz/submissions/hubs_rosestoica.csv

=== Vizualizare completa ===
