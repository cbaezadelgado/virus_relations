#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 12:25:43 2022

@author: Carlos Baeza-Delgado
"""

import os
import sys
import json
import pandas as pd
import yaml


def is_desired_relation(sent: dict,
                        triplet: dict,
                        relations: list):
    '''
    This function checks if a relation between a virus and a gene exists in a
    given sentence.

    Parameters
    ----------
    sent : dict
        A dictionary representing a sentence and its associated information,
        including a list of mentions.
    triplet : dict
        A dictionary representing a relation triplet, including mentions, the
        indices of the head and tail mentions and the relation.

    Returns
    -------
    bool
        True if a relation between a virus and a gene exists in the sentence,
        otherwise False.
    '''

    types = set()
    types.add(sent['mention_list'][triplet['h_mention']]['type'])
    types.add(sent['mention_list'][triplet['t_mention']]['type'])
    desired_relation = False
    for relation in relations:
        relation = set(relation)
        if relation == types:
            desired_relation = True
            break
    return desired_relation


def get_real_mentions(sent: dict,
                      triplet: dict):
    '''
    This function retrieves the real names of a virus and a gene from a
    sentence and a triplet (relation between two different mentions).

    Parameters
    ----------
    sent : dict
        A dictionary representing a sentence and its associated information,
        including a list of mentions.
    triplet : dict
        A dictionary representing a relation triplet, including mentions, the
        indices of the head and tail mentions and the relation.

    Returns
    -------
    tuple
        A tuple of two strings, representing the real names of a virus and a
        gene in the sentence.
    '''

    if sent['mention_list'][triplet['h_mention']]['type'] == 'virus':
        virus = sent['mention_list'][triplet['h_mention']]['name']
        other = sent['mention_list'][triplet['t_mention']]['name']
        virus_sci = sent['mention_list'][triplet['h_mention']]['ScientificName']
        other_type = sent['mention_list'][triplet['t_mention']]['type']
        try:
            other_sci = sent['mention_list'][triplet['t_mention']]['ScientificName']
        except KeyError:
            if sent['mention_list'][triplet['t_mention']]['id'] != ['CUI-less']:
                other_sci = sent['mention_list'][triplet['t_mention']]['id'][0]
            else:
                other_sci = ''
    else:
        virus = sent['mention_list'][triplet['t_mention']]['name']
        other = sent['mention_list'][triplet['h_mention']]['name']
        virus_sci = sent['mention_list'][triplet['t_mention']]['ScientificName']
        other_type = sent['mention_list'][triplet['h_mention']]['type']
        try:
            other_sci = sent['mention_list'][triplet['h_mention']]['ScientificName']
        except KeyError:
            if sent['mention_list'][triplet['h_mention']]['id'] != ['CUI-less']:
                other_sci = sent['mention_list'][triplet['h_mention']]['id'][0]
            else:
                other_sci = ''
    return virus, virus_sci, other, other_sci, other_type


def get_df_relations(path: str,
                     file_name: str,
                     relations: list):
    '''
    This function reads in a JSON file of sentences (given by path and
    file_name parameters), checks if any relation between a virus and a gene
    exists, and returns a dataframe of the relations.

    Parameters
    ----------
    path : str
        The file path to the JSON file.
    file_name : str
        The name of the JSON file.

    Returns
    -------
    df_rel : pandas.DataFrame
        A dataframe of relations between viruses and genes, including the PMID,
        sentence ID, first and second terms of the relation, the relation
        itself, if it is perfect match, and the virus and gene names.
    '''

    with open(os.path.join(path, file_name), 'r', encoding='utf-8') as file:
        sent_list = json.load(file)

    triplets = []

    for sent in sent_list:
        if len(sent['triplet_list']) > 0:
            for tri in sent['triplet_list']:
                if is_desired_relation(sent, tri, relations):
                    (virus_mention,
                     virus_sci,
                     other_mention,
                     other_sci,
                     other_type) = get_real_mentions(sent, tri)
                    virus = f'{virus_mention} ({virus_sci})'
                    if other_sci == '':
                        other = f'{other_mention}'
                    else:
                        other = f'{other_mention} ({other_sci.split(";")[0]})'
                    triplets.append((sent['pmid'],
                                     sent['sent_id'],
                                     tri['triplet'],
                                     virus,
                                     other,
                                     other_type,
                                     sent['sentence']))
                # else:
                #     print(tri['triplet'])

    # unique_triplets = set(['_'.join(x) for x in triplets])

    # unique_triplets = [x.split('_') for x in unique_triplets]

    df_rel = pd.DataFrame()

    df_rel['pmid'] = [x[0] for x in triplets]  # pmid
    df_rel['sent_id'] = [x[1] for x in triplets]  # sent_id
    df_rel['term_1'] = [x[2][0] for x in triplets]  # h_term
    df_rel['relation'] = [x[2][1] for x in triplets]  # rel
    df_rel['term_2'] = [x[2][2] for x in triplets]  # t_term
    df_rel['virus'] = [x[3] for x in triplets]  # virus
    df_rel['other_entity'] = [x[4] for x in triplets]  # gene
    df_rel['relation_with'] = [x[5] for x in triplets]  # other_type
    df_rel['sentence'] = [x[6] for x in triplets]  # sentence

    return df_rel


def get_relations(config, spacy_file, openie_file):

    PATH = config['path']

    relations = config['relations']

    df_spacy = get_df_relations(os.path.join(PATH, 'temp'),
                                spacy_file,
                                relations)
    df_spacy.insert(8, 'source', 'Spacy')

    df_openie = get_df_relations(os.path.join(PATH, 'temp'),
                                 openie_file,
                                 relations)
    df_openie.insert(8, 'source', 'OpenIE')

    for df in [df_spacy, df_openie]:
        df['pmid'] = df['pmid'].astype(int)
        df.sort_values(['pmid', 'sent_id', 'virus', 'other_entity'],
                       inplace=True)

    df_all = pd.concat([df_openie, df_spacy], ignore_index=True)
    df_all['pmid'] = df_all['pmid'].astype(int)
    df_all.sort_values(['pmid', 'sent_id', 'virus', 'other_entity'],
                       inplace=True)

    df_spacy.to_excel(os.path.join(PATH,
                                   'temp',
                                   spacy_file.replace('.json', '.xlsx')),
                      index=False)

    df_openie.to_excel(os.path.join(PATH,
                                    'temp',
                                    openie_file.replace('.json', '.xlsx')),
                       index=False)
    return df_all
