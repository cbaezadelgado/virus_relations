i#!/bin/bash

set -eux

python openie_main.py \
--source_file $1 \
--target_dir ./target_dir \
--jar_dir ./jar_dir/ \
--batch_size 2000 \
--batch_processes 4 \
--process_memory 6g \
--indent 2 \
2>&1 | tee log.txt
