#!/bin/bash

set -eux

python spacy_main.py \
--source_file $1 \
--batch_size 50000 \
--indent 2 \
2>&1 | tee log.txt

