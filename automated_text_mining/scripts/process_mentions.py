#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 12:45:05 2022

@author: Carlos Baeza-Delgado
"""


import pandas as pd
from alive_progress import alive_bar
from scripts.entrez_functions import get_entrez_tax_info
from multiprocessing import Pool


def is_virus_specie(hit):
    virus_specie = False
    if hit:
        if 'Division' in hit:
            if hit['Division'] == 'Viruses':
                if hit['Rank'] == 'species':
                    virus_specie = True
                elif hit['Rank'] == 'no rank' and 'LineageEx' in hit.keys():
                    if hit['LineageEx'][-1]['Rank'] == 'species':
                        virus_specie = True
                elif hit['Rank'] == 'serotype' and 'LineageEx' in hit.keys():
                    for taxon in hit['LineageEx']:
                        if taxon['Rank'] == 'species':
                            virus_specie = True
            elif hit['Division'] == 'Synthetic and Chimeric':
                if 'viruses' in hit['Lineage'].lower():
                    virus_specie = True
    return virus_specie


def get_specie(taxid: str,
               entrez_specie: dict):
    if taxid in entrez_specie:
        return entrez_specie[taxid]
    else:
        return None
        # return get_entrez_tax_info(taxid)[0]


def manage_norm_non_alternative(taxid, annotation, entrez_specie):
    specie = get_specie(taxid, entrez_specie)
    if is_virus_specie(specie):
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = specie['ScientificName']
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['alternative'] = 'None'
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def manage_norm_single_alternative(taxid,
                                   new_taxid,
                                   annotation,
                                   entrez_specie):
    specie = get_specie(taxid, entrez_specie)
    if new_taxid == taxid:
        # print('\t\t- Alternative ID = original ID...')
        if is_virus_specie(specie):
            # print('\t\t\t- IS VIRUS SPECIE!!')
            annotation['ScientificName'] = specie['ScientificName']
            annotation['type'] = 'virus'
            annotation['checked'] = True
            annotation['alternative'] = 'single_equal'
        else:
            # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
            annotation['checked'] = True
    else:
        # print('\t\t\t- NOT SAME ID!')
        new_specie = entrez_specie[new_taxid]
        if specie:
            if is_virus_specie(specie) and is_virus_specie(new_specie):
                name = f'{specie["ScientificName"]};{new_specie["ScientificName"]}'
            elif is_virus_specie(specie) and not is_virus_specie(new_specie):
                name = f'{specie["ScientificName"]}'
            elif not is_virus_specie(specie) and is_virus_specie(new_specie):
                name = f'{new_specie["ScientificName"]}'
            else:
                name = False
        else:
            if is_virus_specie(new_specie):
                name = f'{new_specie["ScientificName"]}'
            else:
                name = False

        if name:
            annotation['ScientificName'] = specie['ScientificName']
            annotation['type'] = 'virus'
            annotation['checked'] = True
            annotation['alternative'] = 'single_diff'
        else:
            # print('\t\t\t\t--> NONE VIRUS SPECIE, REMOVING!!')
            annotation['checked'] = True


def manage_norm_multiple_alternative(id_list,
                                     taxid,
                                     annotation,
                                     entrez_specie):
    specie = get_specie(taxid, entrez_specie)
    original_name = None
    names = set()
    if is_virus_specie(specie):
        original_name = specie['ScientificName']
        names.add(original_name)
    for new_taxid in id_list:
        new_specie = entrez_specie[new_taxid]
        if is_virus_specie(new_specie):
            new_name = new_specie['ScientificName']
            if new_name not in original_name and original_name not in new_name:
                names.add(new_name)
    if None in names:
        names.remove(None)
    if len(names) == 1:
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = list(names)[0]
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['alternative'] = 'single'
    elif len(names) > 1:
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = ';'.join(names)
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['alternative'] = 'multiple'
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def parse_normalized_species(annotation, sp_raw_to_ids, entrez_specie):
    # print(f'......parsing mention: {annotation["mention"]}.......')
    taxid = annotation['id'][0].split(':')[-1]
    alternative = sp_raw_to_ids[annotation['mention']]

    if len(alternative) == 0:
        # print('\t- NO ALTERNATIVE, checing_virus_specie...')
        manage_norm_non_alternative(taxid, annotation, entrez_specie)
    elif len(alternative) == 1:
        new_taxid = alternative[0]
        manage_norm_single_alternative(taxid, new_taxid, annotation,
                                       entrez_specie)
    else:
        # print('\t\t\t-  MORE THAN ONE ID')
        id_list = alternative
        manage_norm_multiple_alternative(id_list, taxid, annotation,
                                         entrez_specie)


def manage_not_norm_not_matched(annotation,
                                alias_virus,
                                sp_raw_to_ids,
                                entrez_specie):
    new_mention = None
    # mention = annotation['mention']
    for virus, aliases in alias_virus.items():
        if annotation['mention'].lower() in aliases:
            # print(f'\tNOT NORM {mention} ALIAS MATCHED')
            new_mention = virus
    if new_mention:
        id_list = sp_raw_to_ids[new_mention]
        if len(id_list) == 0:
            annotation['alias_matched'] = False
        elif len(id_list) == 1:
            # print('\tNOT NORM NOT MATCHED ALIAS MATCHED SINGLE')
            taxid = id_list[0]
            annotation['alias_matched'] = True
            manage_not_norm_single(taxid, annotation, entrez_specie)
        else:
            # print('\tNOT NORM NOT MATCHED ALIAS MATCHED MULTIPLE')
            annotation['alias_matched'] = True
            manage_not_norm_multiple(id_list, annotation, entrez_specie)
    else:
        # print(f'\tNOT NORM {mention} NOT ENTREZ NOT ALIAS MATCHED')
        annotation['checked'] = True


def manage_not_norm_single(taxid,
                           annotation,
                           entrez_specie):
    specie = get_specie(taxid, entrez_specie)
    if is_virus_specie(specie):
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = specie['ScientificName']
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['normalization'] = 'single'
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def manage_not_norm_multiple(id_list,
                             annotation,
                             entrez_specie):
    names = set()
    for new_id in id_list:
        new_specie = entrez_specie[new_id]
        if is_virus_specie(new_specie):
            names.add(new_specie['ScientificName'])

    if len(names) == 1:
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = list(names)[0]
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['normalization'] = 'single'
    if len(names) > 1:
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = ';'.join(names)
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['normalization'] = 'multiple'
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def parse_not_normalized_species(annotation,
                                 alias_virus,
                                 sp_raw_to_ids,
                                 entrez_specie):
    # print(f'......parsing NOT NORM mention: {annotation["mention"]}.......')
    mention = annotation['mention']
    annotation['alias_matched'] = False
    if mention in sp_raw_to_ids:
        id_list = sp_raw_to_ids[mention]
    else:
        id_list = []

    if len(id_list) == 0:
        # print(f'\tNOT NORM {mention} NOT MATCHED')
        manage_not_norm_not_matched(annotation, alias_virus, sp_raw_to_ids,
                                    entrez_specie)
    elif len(id_list) == 1:
        # print('\tNOT NORM MATCHED')
        taxid = id_list[0]
        manage_not_norm_single(taxid, annotation, entrez_specie)
    else:
        manage_not_norm_multiple(id_list, annotation, entrez_specie)


def parse_species(annotation,
                  alias_virus,
                  sp_raw_to_ids,
                  entrez_specie):
    if annotation['id'] != ['CUI-less']:
        parse_normalized_species(annotation, sp_raw_to_ids, entrez_specie)
    else:
        parse_not_normalized_species(annotation, alias_virus, sp_raw_to_ids,
                                     entrez_specie)


def manage_disease_not_matched(annotation,
                               alias_virus,
                               sp_raw_to_ids,
                               entrez_specie):
    # mention = annotation['mention']
    new_mention = None
    for virus, aliases in alias_virus.items():
        if annotation['mention'].lower() in ([virus] + aliases):
            # print(f'\tDISEASE {mention} ALIAS MATCHED')
            new_mention = virus
    if new_mention:
        annotation['alias_matched'] = False
        id_list = sp_raw_to_ids[new_mention]
        if len(id_list) == 0:
            # print(f'\tDISEASE {mention} ALIAS MATCHED NOT ENTREZ')
            annotation['checked'] = True
        elif len(id_list) == 1:
            # print('\tDISEASE NOT MATCHED ALIAS MATCHED SINGLE')
            taxid = id_list[0]
            manage_disease_single(taxid, annotation, entrez_specie)
            annotation['alias_matched'] = True
        else:
            # print('\tDISEASE NOT MATCHED ALIAS MATCHED MULTIPLE')
            manage_disease_multiple(id_list, annotation, entrez_specie)
            annotation['alias_matched'] = True
    else:
        # print(f'\tNOT NORM {mention} NOT ENTREZ NOT ALIAS MATCHED')
        annotation['checked'] = True


def manage_disease_single(taxid: str,
                          annotation: dict,
                          entrez_specie: dict):

    specie = get_specie(taxid, entrez_specie)
    if is_virus_specie(specie):
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = specie['ScientificName']
        annotation['id'] = ['CUI-less']
        annotation['type'] = 'virus'
        annotation['checked'] = True
        annotation['normalization'] = 'disease_single'
        annotation['alias_matched'] = False
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def manage_disease_multiple(id_list,
                            annotation,
                            entrez_specie):
    names = set()
    for new_id in id_list:
        new_search = entrez_specie[new_id]
        if is_virus_specie(new_search):
            names.add(new_search['ScientificName'])

    if len(names) > 0:
        # print('\t\t\t- IS VIRUS SPECIE!!')
        annotation['ScientificName'] = ';'.join(names)
        annotation['type'] = 'virus'
        annotation['id'] = ['CUI-less']
        annotation['checked'] = True
        annotation['normalization'] = 'disease_multiple'
        annotation['alias_matched'] = False
    else:
        # print('\t\t\t- NOT VIRUS SPECIE, REMOVING!!')
        annotation['checked'] = True


def parse_disease(annotation,
                  alias_virus,
                  sp_raw_to_ids,
                  entrez_specie):
    mention = annotation['mention']
    # print(f'......parsing DISEASE mention: {mention}.......')
    if mention in sp_raw_to_ids:
        id_list = sp_raw_to_ids[mention]
    else:
        id_list = []
    if len(id_list) == 0:
        # print('\tNOT MATCHED')
        annotation['checked'] = True
        # manage_disease_not_matched(annotation, alias_virus, sp_raw_to_ids,
        #                            entrez_specie)
    elif len(id_list) == 1:
        # print('\tMATCHED')
        taxid = id_list[0]
        manage_disease_single(taxid, annotation, entrez_specie)
    else:
        # print('\tMATCHED')
        manage_disease_multiple(id_list, annotation, entrez_specie)


def parse_gene(annotation,
               gene_raw_to_ids,
               entrez_genes):
    gene_id = annotation['id'][0].split(':')[-1]
    name = entrez_genes[gene_id]['gene_desc']
    names = [name]
    mention = annotation['mention']
    if mention in gene_raw_to_ids:
        id_alt = gene_raw_to_ids[mention]
    else:
        id_alt = False
    if id_alt and id_alt != gene_id:
        name_alt = entrez_genes[id_alt]['gene_desc']
        names.append(name_alt)
    if len(names) > 0:
        annotation['ScientificName'] = ';'.join(names)
    annotation['checked'] = True


def parse_mentions(args: tuple):
    abstract = args[0]
    alias_virus = args[1]
    sp_raw_to_ids = args[2]
    gene_raw_to_ids = args[3]
    entrez_species = args[4]
    entrez_genes = args[5]
    parse_mentions = args[6]
    for annotation in abstract['mention_list']:
        annotation['name'] = str(annotation['mention'])
        if annotation['type'] == 'species' and 'species' in parse_mentions:
            parse_species(annotation, alias_virus, sp_raw_to_ids,
                          entrez_species)
        elif annotation['type'] == 'disease' and 'disease' in parse_mentions:
            parse_disease(annotation, alias_virus, sp_raw_to_ids,
                          entrez_species)
        elif annotation['type'] == 'gene' and 'genes' in parse_mentions:
            parse_gene(annotation, gene_raw_to_ids, entrez_genes)
        else:
            annotation['name'] = annotation['mention']

    return abstract


def parse_all_mentions(annotated_abstracts,
                       alias_virus,
                       sp_raw_to_ids,
                       gene_raw_to_ids,
                       entrez_species,
                       entrez_genes):

    total_mentions = sum([len(x['mention_list']) for x in annotated_abstracts])

    annotated_parsed = []
    print('\n- Processing mentions of annotated sentences')
    with alive_bar(total_mentions) as bar:
        for abstract in annotated_abstracts[:]:
            new_abstract = {k: v for k, v in abstract.items()
                            if k != 'mention_list'}
            new_mention_list = [{k: v for k, v in ment.items()}
                                for ment in abstract['mention_list']]
            new_abstract['mention_list'] = new_mention_list
            annotated_parsed.append(new_abstract)
        for abstract in annotated_parsed:
            for annotation in abstract['mention_list']:
                bar()
                if annotation['type'] == 'species':
                    parse_species(annotation, alias_virus, sp_raw_to_ids,
                                  entrez_species)
                elif annotation['type'] == 'disease':
                    parse_disease(annotation, alias_virus, sp_raw_to_ids,
                                  entrez_species)
                elif annotation['type'] == 'gene':
                    parse_gene(annotation, gene_raw_to_ids, entrez_genes)
                else:
                    annotation['name'] = annotation['mention']

    return annotated_parsed


def parse_all_mentions_parallel(annotated_abstracts,
                                alias_virus,
                                sp_raw_to_ids,
                                gene_raw_to_ids,
                                entrez_species,
                                entrez_genes,
                                mentions_to_parse,
                                n_batch,
                                total_batches,
                                n_proc=10):

    annotated_parsed = []
    for abstract in annotated_abstracts[:]:
        new_abstract = {k: v for k, v in abstract.items()
                        if k != 'mention_list'}
        if len(abstract['mention_list']) > 0:
            new_mention_list = [{k: v for k, v in m.items()}
                                for m in abstract['mention_list']]
        else:
            new_mention_list = []
        new_abstract['mention_list'] = new_mention_list
        annotated_parsed.append(new_abstract)

    args = [(abstract,
             alias_virus,
             sp_raw_to_ids,
             gene_raw_to_ids,
             entrez_species,
             entrez_genes,
             mentions_to_parse) for abstract in annotated_parsed]

    total_mentions = sum([len(x['mention_list']) for x in annotated_parsed])
    annotated_parsed_2 = []
    print(f'\n- Processing mentions of annotated sentences. BATCH {n_batch}/{total_batches}')
    with alive_bar(total_mentions) as bar, Pool(n_proc) as p:
        for res in p.imap(parse_mentions, args):
            if len(res['mention_list']) > 0:
                bar(len(res['mention_list']))
            annotated_parsed_2.append(res)
    p.close()

    return annotated_parsed_2


def get_and_print_total_mentions(abstract_lst):
    total_mentions = {}

    for abstract in abstract_lst:
        for annotation in abstract['mention_list']:
            mention_class = annotation['type']
            if annotation['type'] not in total_mentions:
                total_mentions[mention_class] = 1
            else:
                total_mentions[mention_class] += 1

    df_counts = pd.DataFrame()

    df_counts['mention_class'] = total_mentions.keys()
    df_counts['counts'] = total_mentions.values()
    df_counts.set_index('mention_class', inplace=True)

    print()
    print('MENTIONS SUMMARY')

    print(df_counts)
    return df_counts


def get_numbers_virus(abstracts_lst):
    checked = 0
    checked_not_norm = 0
    alternatives = []
    normalization = []
    normalization_alias = []
    alias_matched = 0

    for abstract in abstracts_lst:
        for ann in abstract['mention_list'][:]:
            if ann['type'] == 'virus':
                if ann['id'] != ['CUI-less']:
                    checked += ann['checked']
                    alternatives.append(ann['alternative'])
                else:
                    checked_not_norm += ann['checked']
                    alias_matched += ann['alias_matched']
                    if ann['alias_matched']:
                        normalization_alias.append(ann['normalization'])
                    else:
                        normalization.append(ann['normalization'])

    print('\n\n**** NORMALIZED SPECIES ****')
    for i in set(alternatives):
        print(f'{i}: {alternatives.count(i)}')

    print('\n\n**** NOT NORMALIZED SPECIES ****')
    print('\nNO ALIAS:')
    for i in set(normalization):
        print(f'{i}: {normalization.count(i)}')
    print('\nALIAS:')
    for i in set(normalization_alias):
        print(f'{i}: {normalization_alias.count(i)}')


def co_ocurrence(sentence: dict,
                 relations: list):
    '''
    Function to identify if in the mentions of an abstract sentence there are
    co-ocurrence of virus and gene mention types.

    Parameters
    ----------
    sentence : dict
        dict (from json) of a sentence, with info and list of mentions. Each
        mention is also a dict with its info.
    co_ocurrences : list
        list of tuples, where each tuple a the co-ocurrence (a pair of
        mentions type) of interest.

    Returns
    -------
    bool
        True if there are any if the co-ocurrences in the sentence, False
        otherwise.

    '''
    ocurrence = False
    for relation in relations:
        first = False
        second = False
        for mention in sentence['mention_list']:
            if mention['type'] == relation[0]:
                first = True
            elif mention['type'] == relation[1]:
                second = True
        ocurrence = first and second
        if ocurrence:
            break
    return ocurrence


def filter_mentions(annotated_abstracts: dict,
                    relations: list):

    # Make a deepcopy of database
    recvir = []

    keep_mentions = set().union(*relations)

    rels_print = ['-'.join(rel) for rel in relations]
    rels_print = ' or '.join(rels_print)

    print('\n- Filtering mentions, keep sentences with co-ocurring ' +
          f'{rels_print} mentions')

    # Remove mentions other than virus and genes
    with alive_bar(len(annotated_abstracts)) as bar:
        for a in annotated_abstracts:
            new_a = {}
            new_a['mention_list'] = []
            for k, v in a.items():
                if k != 'mention_list':
                    new_a[k] = v
            for ann in a['mention_list']:
                if ann['type'] in keep_mentions:
                    new_a['mention_list'].append(ann)
            recvir.append(new_a)
            bar()

    # Make a deepcopy of recvir database
    co_recvir = []

    # Remove sentences with no virus-gene co-ocurrence
    CO = 0
    for sent in recvir:
        # print(co_ocurrence(a))
        CO += co_ocurrence(sent, relations)
        if co_ocurrence(sent, relations):
            co_recvir.append(sent)

    return co_recvir, CO
