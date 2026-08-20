seq=
sp=
ref=
python=

mkdir -p ${sp}
cd ${sp}

${python} ../01_generate_candidates.py --fasta ${seq} --outdir out

${python} ../02_pair_primers.py --fasta ${seq}  --candidates out/candidates_pass.tsv --out out/pairs_pass.tsv
${python} ../02b_rescore_pairs.py -i out/pairs_pass.tsv
${python} ../03_prepare_blast.py --fasta ${seq} --pairs out/pairs_pass.tsv --top 200 --outdir out/blast

makeblastdb -in ${ref} -dbtype nucl -parse_seqids -out ${sp}_db

blastn -task blastn-short -query out/blast/blast_primers.fa -db ${sp}_db -dust no -soft_masking false -evalue 1000 -max_target_seqs 10000 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qseq sseq' -out out/blast/primer_hits.tsv

${python} ../04_parse_blast.py --pairs out/blast/top_pairs.tsv --blast out/blast/primer_hits.tsv --out out/blast/pairs_with_blast.tsv

