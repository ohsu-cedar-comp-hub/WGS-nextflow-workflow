#!/bin/bash

#SBATCH --job-name=download
#SBATCH --partition=exacloud
#SBATCH --account=cedar
#SBATCH --qos=normal
#SBATCH --time=12:00:00
#SBATCH --mem=2GB
#SBATCH --mincpus=1
#SBATCH --cpus-per-task=1
#SBATCH --output=download_%j.out
#SBATCH --error=download_%j.err

# Temporary directory path
TMPDIR="path/to/temp/directory/"
# Bucket path
BUCKET="rgw/coh"
# File containing list of files to download
FILELIST="/path/to/list/of/files.txt"

## copy from bucket
while IFS= read -r FILE
do
  mc cp "${BUCKET}/${FILE}" "${TMPDIR}/${FILE}"
  if [ $? -eq 0 ]; then
    echo "Success: $FILE"
  else
    echo "$FILE failed"
  fi
done < "$FILELIST"


