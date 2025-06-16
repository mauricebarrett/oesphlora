#Set up GTDB classifier

wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/genomic_files_all/ssu_all_r226.fna.gz
unzip ssu_all_r226.fna.gz

grep -E '^>' ssu_all_r226.fna | tr -d '>' | sed 's/ d__/\td__/' | sed 's/\[.*//' >ssu_all_r226_taxonomy.tsv

qiime tools import \
    --input-path ssu_all_r226_taxonomy.tsv \
    --input-format HeaderlessTSVTaxonomyFormat \
    --type 'FeatureData[Taxonomy]' \
    --output-path ssu_all_r226_taxonomy.qza

qiime tools import \
    --input-path ssu_all_r226.fna \
    --input-format DNAFASTAFormat \
    --type 'FeatureData[Sequence]' \
    --output-path ssu_all_r226_sequences.qza

qiime feature-classifier extract-reads \
    --i-sequences ssu_all_r226_sequences.qza \
    --p-f-primer CCTACGGGNGGCWGCAG \
    --p-r-primer GACTACHVGGGTATCTAATCC \
    --p-min-length 200 \
    --p-max-length 600 \
    --p-n-jobs 18 \
    --p-read-orientation 'both' \
    --o-reads ssu_all_r226_sequences_extracted.qza

qiime rescript dereplicate \
    --i-sequences ssu_all_r226_sequences_extracted.qza \
    --i-taxa ssu_all_r226_taxonomy.qza \
    --p-mode 'majority' \
    --o-dereplicated-sequences ssu_all_r226_sequences_extracted_derep.qza \
    --o-dereplicated-taxa ssu_all_r226_taxonomy_extracted_derep.qza

qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads ssu_all_r226_sequences_extracted_derep.qza \
    --i-reference-taxonomy ssu_all_r226_taxonomy_extracted_derep.qza \
    --o-classifier ssu_all_r226_extracted_derep_classifier.qza
