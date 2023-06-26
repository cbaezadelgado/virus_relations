#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 15 14:51:06 2023

@author: Carlos Baeza-Delgado
"""

import os
import sys
import json

PATH = '/home/qhr/virus_relations/databases/'
with open(os.path.join(PATH, 'annotated_abstracts.json'), mode='r') as f:
    stored_annotated_abstracts = json.load(f)

print(f'\n{len(stored_annotated_abstracts)} abstracts annotated')

annotated = []
for pmid, sents in stored_annotated_abstracts.items():
    annotated.append((pmid, sents))

removed = 0
for abstract in annotated:
    mentions = False
    for sent in abstract[1]:
        if len(sent['mention_list']) > 0:
            mentions = True
    if not mentions:
        del stored_annotated_abstracts[abstract[0]]
        removed += 1


with open(os.path.join(PATH, 'annotated_abstracts.json'), mode='w') as f:
    json.dump(stored_annotated_abstracts, f)

print(f'\n{removed} abstracts removed. ({len(stored_annotated_abstracts)} abstracts remaining)')
