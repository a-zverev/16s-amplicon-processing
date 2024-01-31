#!/bin/bash

QIIME_ENV_PATH="/home/alexey/.conda/envs/qiime2/bin/"
PATH="$QIIME_ENV_PATH:$PATH"

FOLDER="tree_tmp"
FASTA_FILE=$1
TREE_FILE="${2:-$'tree.nwk'}"

[ -f "$FASTA_FILE" ] || { echo "$FASTA_FILE NOT FOUND" ; exit 1 ;}

{ [[ -d "$FOLDER" ]] && rm -rf "$FOLDER"; } ; mkdir -p "$FOLDER";
cd $FOLDER

# import references
echo "Import data..."

qiime tools import \
  --input-path ../$FASTA_FILE \
  --output-path rep-seqs.qza \
  --type 'FeatureData[Sequence]'

# alignment and tree construction
echo "Alignment..."
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# export tree
echo "Exporting tree.."
 qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-tree

# move and cleanup
mv ./exported-tree/tree.nwk ./../$TREE_FILE
rm -r ../$FOLDER