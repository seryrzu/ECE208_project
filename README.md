# ECE208_project

To extract:
```
python genes_extractor.py -g ../data/reference_locus/igh_locus_Mar6_18.fa -c ../data/exact_pos/IGHV_exact_pos.csv
```


Stage 1:
```
clustalo -i data/germline_new/Homo_sapiens/ig/IGHV_with_RSS.fa -o data/germline_new/Homo_sapiens/ig/IGHV_with_RSS.aln
raxmlHPC-AVX  -p 21319301 -s data/germline_new/Homo_sapiens/ig/IGHV_with_RSS.aln -n IGHV_tree -m GTRCAT
```
