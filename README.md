# raw_gene_prediction

Script for raw (only blast based) protein prediction.
1. run blast with recommended  blast options: `-outfmt 5 -gapopen 11 -gapextend 2
-dust no -soft_masking`

2. run raw_gene_prediction

## usage
`./raw_prediction.py query_seqs genome blast_output short_name output`

optional arguments:
  --genetic-code GENETIC_CODE
                        altrnative genetic code in a form: "TAA:Q;TAG:Q"
  --threads THREADS     number of threads
  --hits HITS           number of hits for one protein

