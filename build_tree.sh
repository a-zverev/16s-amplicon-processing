#!/bin/bash

set -euo pipefail

QIIME_ENV_PATH="/home/alexey/.conda/envs/qiime2/bin/"
export PATH="$QIIME_ENV_PATH:$PATH"

# Input arguments
FASTA_FILE="$1"
TREE_FILE="${2:-tree.nwk}"

# Temp working folder
FOLDER="tree_tmp"

# Check input
if [[ ! -f "$FASTA_FILE" ]]; then
  echo "ERROR: Input FASTA file '$FASTA_FILE' not found."
  exit 1
fi

# Check qiime is available
if ! command -v qiime &> /dev/null; then
  echo "ERROR: 'qiime' command not found. Check your QIIME2 environment."
  exit 1
fi

# Prepare temporary folder
rm -rf "$FOLDER"
mkdir -p "$FOLDER"
cd "$FOLDER"

# Import sequences
echo "Step 1: Importing sequences..."
qiime tools import \
  --input-path "../$FASTA_FILE" \
  --output-path rep-seqs.qza \
  --type 'FeatureData[Sequence]'

# Align and build tree
echo "Step 2: Building phylogenetic tree..."
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

# Export rooted tree
echo "Step 3: Exporting tree..."
qiime tools export \
  --input-path rooted-tree.qza \
  --output-path exported-tree

# Move tree file
mv exported-tree/tree.nwk "../$TREE_FILE"

# Cleanup
cd ..
rm -rf "$FOLDER"

echo "Tree written to '$TREE_FILE'"
