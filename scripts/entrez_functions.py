#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 11:55:34 2022

@author: Carlos Baeza-Delgado
"""

import re
import os
import json
import time
import http.client
from urllib.error import HTTPError
from Bio import Entrez
from unidecode import unidecode
from multiprocessing import Pool
from alive_progress import alive_bar


def load_json_database(path: str,
                       filename: str):
    '''
    DESCRIPTION.

    Parameters
    ----------
    path : str
        DESCRIPTION.
    filename : str
        DESCRIPTION.

    Returns
    -------
    my_json_dict : dict
        DESCRIPTION.

    '''
    try:
        with open(os.path.join(path, filename),
                  mode='r',
                  encoding='utf-8') as json_file:
            my_json_dict = json.load(json_file)
            n_entries = len(my_json_dict)
            # print(f'\n\nDatabase {filename} loaded successfully from {path}' +
            #       f' with {n_entries} entries')
    except FileNotFoundError:
        my_json_dict = {}  # First run, then load saved json files
        # print(f'\n\nDatabase {filename} NOT found in {path}, ' +
        #      'initializing database with 0 entries')
    return my_json_dict


def update_json_database(path: str,
                         filename: str,
                         database: dict):
    '''
    DESCRIPTION.

    Parameters
    ----------
    path : str
        DESCRIPTION.
    filename : str
        DESCRIPTION.
    database : dict
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    with open(os.path.join(path, filename),
              mode='w',
              encoding='utf-8') as json_file:
        json.dump(database, json_file)
        # n_entries = len(database)
        # print(f'\n\nDatabase {filename} with {n_entries} entries updated ' +
        #       f'successfully in {path}')


def search_entrez_gene_mention(mention, entrez_gene_mentions):
    if mention not in entrez_gene_mentions.keys():
        retrieved = False
        print(f'--> Retrieving info for GENE mention: {mention}',
              end=' ....... ')
        while not retrieved:
            try:
                handle = Entrez.esearch(db='gene',
                                        term=mention,
                                        retmode='xml')
                rec = Entrez.read(handle)
                handle.close()
                print('DONE')
                retrieved = True
            except HTTPError:
                print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
                time.sleep(2)
            except ValueError:
                print('\n****** WAITING FOR ENTREZ, ValueError *******',
                      end='\t')
                time.sleep(2)
            except http.client.IncompleteRead:
                print('\n****** WAITING FOR ENTREZ, IncompleteRead *******',
                      end='\t')
                time.sleep(2)
        entrez_gene_mentions[mention] = rec
    return entrez_gene_mentions[mention]


def extract_paper_info(paper):
    paper_dict = {}

    if len(paper['PubmedArticle']) > 0:
        try:
            pmid = paper['PubmedArticle'][0]['MedlineCitation']['PMID']
        except KeyError:
            pmid = 'None'

        try:
            article = paper['PubmedArticle'][0]['MedlineCitation']['Article']
            title = article['ArticleTitle']
        except KeyError:
            title = 'None'

        try:
            abstract = article['Abstract']['AbstractText'][0]
            abstract = unidecode(re.sub('<[^<]+?>', '', abstract))
        except KeyError:
            abstract = 'None'

        try:
            year = article['Journal']['JournalIssue']['PubDate']['Year']
        except KeyError:
            year = 'None'

        try:
            other_ids = paper['PubmedData']['ArticleIdList']
            pmc = 'None'
            for i in other_ids:
                if i.attributes['IdType'] == 'pmc':
                    pmc = str(i)
        except KeyError:
            pmc = 'None'

        try:
            descriptors = []
            medline_citation = paper['PubmedArticle'][0]['MedlineCitation']
            descriptors_entrez = medline_citation['MeshHeadingList']
            for desc in descriptors_entrez:
                descriptors.append(str(desc['DescriptorName']))
        except KeyError:
            descriptors = 'None'

        paper_dict['pmid'] = str(pmid)
        paper_dict['title'] = str(title)
        paper_dict['abstract'] = str(abstract)
        paper_dict['year'] = str(year)
        paper_dict['pmc'] = str(pmc)
        paper_dict['descriptors'] = descriptors

    return paper_dict


def get_entrez_paper_info(pmid, entrez_abstracts):
    if pmid not in entrez_abstracts.keys():
        retrieved = False
        print(f'--> Retrieving info for PMID: {pmid}', end=' ....... ')
        while not retrieved:
            try:
                handle = Entrez.efetch(db='pubmed',
                                       id=pmid,
                                       retmode='xml')
                rec = Entrez.read(handle)
                handle.close()
                print('DONE')
                retrieved = True
            except HTTPError:
                print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
                time.sleep(2)
            except ValueError:
                print('\n****** WAITING FOR ENTREZ, ValueError *******',
                      end='\t')
                time.sleep(2)
            except http.client.IncompleteRead:
                print('\n****** WAITING FOR ENTREZ, IncompleteRead *******',
                      end='\t')
                time.sleep(2)
        entrez_abstracts[pmid] = extract_paper_info(rec)
    return entrez_abstracts[pmid]


def extract_papers_info(paper):
    paper_dict = {}

    try:
        pmid = paper['MedlineCitation']['PMID']
    except KeyError:
        pmid = 'None'

    try:
        article = paper['MedlineCitation']['Article']
        title = article['ArticleTitle']
    except KeyError:
        title = 'None'

    try:
        abstract = article['Abstract']['AbstractText'][0]
        abstract = unidecode(re.sub('<[^<]+?>', '', abstract))
    except KeyError:
        abstract = 'None'

    try:
        year = article['Journal']['JournalIssue']['PubDate']['Year']
    except KeyError:
        year = 'None'

    try:
        other_ids = paper['PubmedData']['ArticleIdList']
        pmc = 'None'
        for i in other_ids:
            if i.attributes['IdType'] == 'pmc':
                pmc = str(i)
    except KeyError:
        pmc = 'None'

    try:
        descriptors = []
        medline_citation = paper['MedlineCitation']
        descriptors_entrez = medline_citation['MeshHeadingList']
        for desc in descriptors_entrez:
            descriptors.append(str(desc['DescriptorName']))
    except KeyError:
        descriptors = 'None'

    paper_dict['pmid'] = str(pmid)
    paper_dict['title'] = str(title)
    paper_dict['abstract'] = str(abstract)
    paper_dict['year'] = str(year)
    paper_dict['pmc'] = str(pmc)
    paper_dict['descriptors'] = descriptors

    return paper_dict


def get_entrez_papers_info(pmids):
    retrieved = False
    # print(f'--> Retrieving info for {len(pmids)} PMIDS', end=' ....... ')
    while not retrieved:
        try:
            handle = Entrez.efetch(db='pubmed',
                                   id=pmids,
                                   retmode='xml')
            rec = Entrez.read(handle)
            handle.close()
            # print('DONE')
            retrieved = True
        except HTTPError:
            print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            print('\n****** WAITING FOR ENTREZ, ValueError *******',
                  end='\t')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            print('\n****** WAITING FOR ENTREZ, IncompleteRead *******',
                  end='\t')
            time.sleep(0.5)
    papers = []
    for paper in rec['PubmedArticle']:
        paper_info = extract_papers_info(paper)
        papers.append(paper_info)

    return papers


def search_entrez_tax_mention(mention, entrez_mentions):
    if mention not in entrez_mentions.keys():
        retrieved = False
        # print(f'--> Retrieving info for TAX mention: {mention}')
        while not retrieved:
            try:
                handle = Entrez.esearch(db='taxonomy',
                                        term=mention,
                                        retmode='xml')
                rec = Entrez.read(handle)
                handle.close()
                retrieved = True
                # print('DONE')
            except HTTPError:
                # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
                time.sleep(2)
            except ValueError:
                # print('\n****** WAITING FOR ENTREZ, ValueError *******')
                time.sleep(2)
            except http.client.IncompleteRead:
                # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
                time.sleep(2)
        entrez_mentions[mention] = rec
    return entrez_mentions[mention]


def get_batches_pmids(pmids, n_processes):
    batch_size = min(int(len(pmids)/n_processes)+1, 9999)
    cuts = list(range(0, len(pmids), batch_size))

    batches = []

    for i in range(len(cuts)):
        if i == len(cuts)-1:
            subset = pmids[cuts[i]:]
        else:
            subset = pmids[cuts[i]:cuts[i+1]]
        batches.append(subset)

    return batches


def download_abstracts(pmids, n_processes=10):
    print('\n- Downloading abstracts from NCBI Pubmed Database')
    batches = get_batches_pmids(pmids, n_processes)
    abstracts_batches = []
    with alive_bar(len(pmids)) as bar, Pool(n_processes) as p:
        for i in p.imap(get_entrez_papers_info, batches):
            # print(i)
            abstracts_batches.append(i)
            bar(len(i))
    p.close()

    # print('\tDONE')
    abstracts_lst = []

    for batch in abstracts_batches:
        abstracts_lst.extend(batch)

    return abstracts_lst


def extract_annotations_and_ids(annotated_abstracts: list):
    '''
    DESCRIPTION.

    Parameters
    ----------
    annotated_abstracts : list
        DESCRIPTION.

    Returns
    -------
    sp_ids : set
        DESCRIPTION.
    gene_ids : set
        DESCRIPTION.
    sp_mentions : list
        DESCRIPTION.
    gene_mentions : list
        DESCRIPTION.

    '''
    print('\n- Searching mentions in annotated abstracts')
    sp_mentions = set()
    sp_ids = set()
    gene_mentions = set()
    gene_ids = set()

    total_mentions = sum([len(x['mention_list']) for x in annotated_abstracts])

    with alive_bar(total_mentions) as bar:
        for sent in annotated_abstracts:
            for mention in sent['mention_list']:
                bar()
                if mention['type'] == 'species':
                    sp_mentions.add(mention['mention'])
                    if mention['id'] != ['CUI-less']:
                        sp_id = mention['id'][0].split(':')[1]
                        sp_ids.add(sp_id)
                elif mention['type'] == 'disease':
                    sp_mentions.add(mention['mention'])
                elif mention['type'] == 'gene':
                    gene_mentions.add(mention['mention'])
                    if mention['id'] != ['CUI-less']:
                        gene_id = mention['id'][0].split(':')[1]
                        gene_ids.add(gene_id)

    # sp_mentions = list(sp_mentions)
    # gene_mentions = list(gene_mentions)
    # sp_ids = list(sp_ids)
    # gene_ids = list(gene_ids)

    return sp_ids, gene_ids, sp_mentions, gene_mentions


# funciones para procesar menciones para la búsqueda en NCBI databases

def fix_sp_mention(x: str):
    '''
    Cleans and normalizes the mention of a specie (text).
    Performs a series of transformations to remove specific characters,
    replace words and symbols, and remove empty spaces.

    Parameters
    ----------
    x : str
        String that contains the specie mention to be cleaned and normalized.

    Returns
    -------
    x : str or None
        Cleaned and normalized string if the length is greater than 2
        alphanumeric characters.
        Otherwise returns None.
    '''

    brackets_pattern = r" \([a-zA-Z0-9-+ /]{1,}\)|\([^a-zA-Z0-9]{1,}\)"
    x = x.replace('viruses', 'virus')
    x = x.replace('_', '-').replace('|', ' ')
    x = re.sub(brackets_pattern, ' ', x)
    x = x.replace('+', ' ')
    x = re.sub(r'[()\"]', ' ', x)
    # Esto sirve para eliminar el último caracter si no está entre los
    # indicados en la expresión [^a-zA-Z0-9)]
    x = re.sub(r"(?P<k>.)[^a-zA-Z0-9)]{1,}$", r"\g<k>", x)
    x = x.replace(' -', ' ')
    x = re.sub(r" +", ' ', x)
    x = x.strip()

    return x if len(re.sub(r"[^a-zA-Z0-9]", '', x)) > 2 else None


def fix_gene_mention(x: str):
    '''
    Cleans and normalizes the mention of a gene (text).
    Performs a series of transformations to remove specific characters,
    replace words and symbols, and remove empty spaces.

    Parameters
    ----------
    x : str
        String that contains the gene mention to be cleaned and normalized.

    Returns
    -------
    x : str or None
        Cleaned and normalized string of text if the length is greater than 2
        alphanumeric characters.
        Otherwise returns None.
    '''

    brackets_pattern = r" \([a-zA-Z0-9-+ /]{1,}\)|\([^a-zA-Z0-9]{1,}\)"
    x = x.replace('viruses', 'virus')
    x = re.sub(r'[Hh]uman[- ]', '', x)
    x = x.replace('_', '-').replace('|', ' ')
    x = re.sub(brackets_pattern, ' ', x)
    x = x.replace('+', ' ')
    x = re.sub(r'[()\"\']', ' ', x)
    # Esto sirve para eliminar el último caracter si no está entre los
    # indicados en la expresión [a-zA-Z0-9)]
    x = re.sub(r"(?P<k>.)[^a-zA-Z0-9)]{1,}$", r"\g<k>", x)
    x = x.replace(' -', ' ')
    x = re.sub(r" +", ' ', x)
    x = x.strip()

    return x if len(re.sub(r"[^a-zA-Z0-9]", '', x)) > 2 else None


# funciones para obtener diccionario fixed_mention: [raw_mentions list]

def get_raw_and_fixed_sp_mentions(mentions_list: list):
    '''
    DESCRIPTION OF THE FUNCTION

    Parameters
    ----------
    mentions_list : list
        DESCRIPTION.

    Returns
    -------
    fixed_to_raw : dict
        DESCRIPTION.

    '''
    fixed_to_raw = dict()
    for raw in mentions_list:
        fixed = fix_sp_mention(raw)
        if fixed:
            if fixed in fixed_to_raw:
                fixed_to_raw[fixed].append(raw)
            else:
                fixed_to_raw[fixed] = [raw]
    return fixed_to_raw


def get_raw_and_fixed_gene_mentions(mentions_list: list):
    '''
    DESCRIPTION OF THE FUNCTION

    Parameters
    ----------
    mentions_list : list
        DESCRIPTION.

    Returns
    -------
    fixed_to_raw : dict
        DESCRIPTION.

    '''
    fixed_to_raw = dict()
    for raw in mentions_list:
        fixed = fix_gene_mention(raw)
        if fixed:
            if fixed in fixed_to_raw:
                fixed_to_raw[fixed].append(raw)
            else:
                fixed_to_raw[fixed] = [raw]
    return fixed_to_raw


def get_ids_from_sp_mention(mention: str):
    '''
    DESCRIPTION.

    Parameters
    ----------
    mention : str
        DESCRIPTION.

    Returns
    -------
    res : list of tuples
        DESCRIPTION.

    '''

    retrieved = False
    while not retrieved:
        try:
            handle = Entrez.esearch(db='taxonomy',
                                    term=mention,
                                    retmode='xml')
            rec = Entrez.read(handle)
            handle.close()
            retrieved = True
        except HTTPError:
            # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            # print('\n****** WAITING FOR ENTREZ, ValueError *******')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
            time.sleep(0.5)
    res = (mention, [str(x) for x in rec['IdList']])

    return res


def get_ids_from_gene_mention(mention: str):
    '''
    DESCRIPTION.

    Parameters
    ----------
    mention : str
        DESCRIPTION.

    Returns
    -------
    res : list of tuples
        DESCRIPTION.

    '''

    retrieved = False
    term = f'(((({mention}[All Fields]) AND human[Organism]) AND "genetype '
    term += 'protein coding"[Properties]) AND "current only"[Filter])'
    while not retrieved:
        try:
            handle = Entrez.esearch(db='gene',
                                    term=term,
                                    retmode='xml',
                                    retmax=1000,
                                    sort='relevance')
            rec = Entrez.read(handle)
            handle.close()
            retrieved = True
        except HTTPError:
            # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            # print('\n****** WAITING FOR ENTREZ, ValueError *******')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
            time.sleep(0.5)
    res = (mention, [str(x) for x in rec['IdList']])

    return res


def download_ids_from_mentions(fixed_to_raw: dict,
                               db: str,
                               n_proc: int = 5):
    '''
    DESCRIPTION.

    Parameters
    ----------
    fixed_to_raw : dict
        DESCRIPTION.
    n_proc : int, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    raw_to_ids : dict
        DESCRIPTION.

    '''

    valid_db = {'species', 'gene'}
    if db not in valid_db:
        raise ValueError(f'db parameter must be one of {valid_db}.')
    if db == 'species':
        get_func = get_ids_from_sp_mention
        db_print = 'Taxonomy'
    else:
        get_func = get_ids_from_gene_mention
        db_print = 'Gene'

    mentions = fixed_to_raw.keys()
    results = []

    print(f'\n- Searching IDs for {db} mentions in NCBI {db_print} Database:')
    with alive_bar(len(mentions)) as bar, Pool(n_proc) as p:
        for i in p.imap(get_func, mentions):
            # print(i)
            results.append(i)
            bar()
    p.close()

    raw_to_ids = dict()
    for result in results:
        fixed, ids = result
        for raw in fixed_to_raw[fixed]:

            if db == 'species':
                raw_to_ids[raw] = ids
            else:
                if len(ids) > 0:
                    raw_to_ids[raw] = ids[0]
                else:
                    raw_to_ids[raw] = ids

    return raw_to_ids


def extract_sp_and_gene_ids(annotated_abstracts: list,
                            alias_virus: dict):
    '''
    DESCRIPTION.

    Parameters
    ----------
    annotated_abstracts : list
        DESCRIPTION.

    Returns
    -------
    sp_ids : set
        DESCRIPTION.
    sp_raw_to_ids : dict
        DESCRIPTION.
    gene_ids : set
        DESCRIPTION.
    gene_raw_to_ids : dict
        DESCRIPTION.

    '''

    # extract mentions and ids from annotated abstracts
    (sp_ids,
     gene_ids,
     sp_mentions,
     gene_mentions) = extract_annotations_and_ids(annotated_abstracts)

    # fix mentions
    sp_fixed_to_raw = get_raw_and_fixed_sp_mentions(sp_mentions)
    gene_fixed_to_raw = get_raw_and_fixed_gene_mentions(gene_mentions)

    # download ids from mentions
    sp_raw_to_ids = download_ids_from_mentions(sp_fixed_to_raw, 'species')
    gene_raw_to_ids = download_ids_from_mentions(gene_fixed_to_raw, 'gene')

    # combine mentions
    sp_ids_mentions = {x for lst in sp_raw_to_ids.values() for x in lst if len(lst) > 0}
    sp_ids = sp_ids.union(sp_ids_mentions)
    gene_ids_mentions = {x for x in gene_raw_to_ids.values() if x != []}
    gene_ids = gene_ids.union(gene_ids_mentions)

    sp_alias_to_raw = dict()
    for raw, ids in sp_raw_to_ids.items():
        if len(ids) == 0:
            for virus, aliases in alias_virus.items():
                if raw.lower() in ([virus] + aliases):
                    if virus in sp_alias_to_raw:
                        sp_alias_to_raw[virus].append(raw)
                    else:
                        sp_alias_to_raw[virus] = [raw]

    sp_raw_to_ids_alias = download_ids_from_mentions(sp_alias_to_raw,
                                                     'species')

    for raw, ids in sp_raw_to_ids_alias.items():
        sp_raw_to_ids[raw] = ids

    return sp_ids, sp_raw_to_ids, gene_ids, gene_raw_to_ids


def extract_sp_ids(sp_ids: set,
                   sp_mentions: set,
                   alias_virus: dict):
    '''
    DESCRIPTION.

    Parameters
    ----------
    annotated_abstracts : list
        DESCRIPTION.

    Returns
    -------
    sp_ids : set
        DESCRIPTION.
    sp_raw_to_ids : dict
        DESCRIPTION.
    gene_ids : set
        DESCRIPTION.
    gene_raw_to_ids : dict
        DESCRIPTION.

    '''

    # fix mentions
    sp_fixed_to_raw = get_raw_and_fixed_sp_mentions(sp_mentions)

    # download ids from mentions
    sp_raw_to_ids = download_ids_from_mentions(sp_fixed_to_raw, 'species')

    # combine mentions
    sp_ids_mentions = {x for lst in sp_raw_to_ids.values() for x in lst if len(lst) > 0}
    sp_ids = sp_ids.union(sp_ids_mentions)

    sp_alias_to_raw = dict()
    for raw, ids in sp_raw_to_ids.items():
        if len(ids) == 0:
            for virus, aliases in alias_virus.items():
                if raw.lower() in ([virus] + aliases):
                    if virus in sp_alias_to_raw:
                        sp_alias_to_raw[virus].append(raw)
                    else:
                        sp_alias_to_raw[virus] = [raw]

    sp_raw_to_ids_alias = download_ids_from_mentions(sp_alias_to_raw,
                                                     'species')

    for raw, ids in sp_raw_to_ids_alias.items():
        sp_raw_to_ids[raw] = ids

    return sp_ids, sp_raw_to_ids


def extract_gene_ids(gene_ids: set,
                     gene_mentions: set):
    '''

    Parameters
    ----------
    gene_ids : set
        DESCRIPTION.
    gene_mentions : set
        DESCRIPTION.

    Returns
    -------
    sp_ids : TYPE
        DESCRIPTION.
    sp_raw_to_ids : TYPE
        DESCRIPTION.

    '''

    # fix mentions
    gene_fixed_to_raw = get_raw_and_fixed_gene_mentions(gene_mentions)

    # download ids from mentions
    gene_raw_to_ids = download_ids_from_mentions(gene_fixed_to_raw, 'gene')

    # combine mentions
    gene_ids_mentions = {x for x in gene_raw_to_ids.values() if x != []}
    gene_ids = gene_ids.union(gene_ids_mentions)

    return gene_ids, gene_raw_to_ids


def extract_gene_info(gene: dict):
    gene_parsed = {}

    gene_id = gene['Entrezgene_track-info']['Gene-track']['Gene-track_geneid']

    try:
        bio_org = gene['Entrezgene_source']['BioSource']['BioSource_org']
        org_ref = bio_org['Org-ref']
        db_tag = org_ref['Org-ref_db'][0]['Dbtag_tag']
        tax_id = db_tag['Object-id']['Object-id_id']
    except KeyError:
        tax_id = 'None'
    try:
        tax_name = org_ref['Org-ref_taxname']
    except KeyError:
        tax_name = 'None'

    try:
        gene_desc = gene['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
    except KeyError:
        gene_desc = 'None'
    try:
        gene_syns = gene['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
    except KeyError:
        gene_syns = 'None'

    try:
        prot_desc = gene['Entrezgene_prot']['Prot-ref']['Prot-ref_desc']
    except KeyError:
        prot_desc = 'None'
    try:
        prot_names = gene['Entrezgene_prot']['Prot-ref']['Prot-ref_name']
    except KeyError:
        prot_names = 'None'

    gene_parsed['gene_id'] = str(gene_id)
    gene_parsed['gene_desc'] = str(gene_desc)
    gene_parsed['gene_syn'] = [str(x) for x in gene_syns]
    gene_parsed['prot_desc'] = str(prot_desc)
    gene_parsed['prot_syn'] = [str(x) for x in prot_names]
    gene_parsed['tax_name'] = str(tax_name)
    gene_parsed['tax_id'] = str(tax_id)

    return gene_parsed


def get_entrez_genes_info(gene_ids: list):
    '''
    DESCRIPTION.

    Parameters
    ----------
    gene_ids : list
        DESCRIPTION.

    Returns
    -------
    genes : list
        DESCRIPTION.

    '''
    retrieved = False
    # print(f'--> Retrieving info for batch size {len(gene_ids)}',
    # end=' ....... ')
    while not retrieved:
        try:
            handle = Entrez.efetch(db='gene',
                                   id=gene_ids,
                                   retmode='xml')
            rec = Entrez.read(handle)
            handle.close()
            # print('DONE')
            retrieved = True
        except HTTPError:
            # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            # print('\n****** WAITING FOR ENTREZ, ValueError *******')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
            time.sleep(0.5)

    genes = []
    for gene in rec:
        gene_info = extract_gene_info(gene)
        genes.append(gene_info)

    return genes


def get_entrez_gene_info(gene_id: str):
    '''
    DESCRIPTION.

    Parameters
    ----------
    gene_ids : list
        DESCRIPTION.

    Returns
    -------
    genes : list
        DESCRIPTION.

    '''
    retrieved = False
    # print(f'--> Retrieving info for GENE ID: {gene_id}', end=' ....... ')
    while not retrieved:
        try:
            handle = Entrez.efetch(db='gene',
                                   id=gene_id,
                                   retmode='xml')
            rec = Entrez.read(handle)
            handle.close()
            # print('DONE')
            retrieved = True
        except HTTPError:
            # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            # print('\n****** WAITING FOR ENTREZ, ValueError *******')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
            time.sleep(0.5)

    gene_info = extract_gene_info(rec[0])

    return gene_info


def get_batches_ids(ids: list,
                    size: int = None,
                    max_size: int = 100,
                    n_processes: int = 5):
    '''
    DESCRIPTION.

    Parameters
    ----------
    ids : list
        DESCRIPTION.
    size : int, optional
        DESCRIPTION. The default is None.
    max_size : int, optional
        DESCRIPTION. The default is 100.
    n_processes : int, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    batches : list of lists
        DESCRIPTION.

    '''
    if size:
        batch_size = size
    else:
        batch_size = min(int(len(ids)/n_processes)+1, max_size)
    cuts = list(range(0, len(ids), batch_size))

    batches = []

    for i in range(len(cuts)):
        if i == len(cuts)-1:
            subset = ids[cuts[i]:]
        else:
            subset = ids[cuts[i]:cuts[i+1]]
        batches.append(subset)

    return batches


def download_genes(ids: list,
                   size: int = None,
                   max_size: int = 50,
                   n_processes: int = 5):
    '''
    DESCRIPTION.

    Parameters
    ----------
    ids : list
        DESCRIPTION.
    size : int, optional
        DESCRIPTION. The default is None.
    max_size : int, optional
        DESCRIPTION. The default is 50.
    n_processes : int, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    genes_info_dict : dict
        DESCRIPTION.

    '''
    print('\n- Downloading genes info from NCBI Gene Database')
    batches = get_batches_ids(ids,
                              size=size,
                              max_size=max_size,
                              n_processes=n_processes)
    genes_info = []
    with alive_bar(len(ids)) as bar, Pool(n_processes) as p:
        for i in p.imap(get_entrez_genes_info, batches):
            # print(i)
            genes_info.append(i)
            bar(len(i))
    p.close()

    # print('\tDONE')
    genes_info_dict = dict()

    for batch in genes_info:
        for gene in batch:
            genes_info_dict[gene['gene_id']] = gene

    return genes_info_dict


def extract_tax_info(tax: dict):
    '''
    DESCRIPTION

    Parameters
    ----------
    tax : dict
        DESCRIPTION.

    Returns
    -------
    tax_parsed : dict
        DESCRIPTION.

    '''
    tax_parsed = {}

    fields_to_keep = ['TaxId',
                      'ScientificName',
                      'Rank',
                      'Division',
                      'Lineage',
                      'LineageEx',
                      'AkaTaxIds'
                      ]

    for field in fields_to_keep:
        try:
            if field == 'LineageEx':
                lin = tax[field]
                new_lin = [{str(k): str(v) for k, v in x.items()} for x in lin]
                tax_parsed[field] = new_lin
            elif field == 'AkaTaxIds':
                tax_parsed[field] = ';'.join(tax['AkaTaxIds'])
            else:
                tax_parsed[field] = str(tax[field])
        except KeyError:
            tax_parsed[field] = None

    return tax_parsed


def get_entrez_tax_info(taxids: list):
    '''
    DESCRIPTION.

    Parameters
    ----------
    taxid : list
        DESCRIPTION.

    Returns
    -------
    species : TYPE
        DESCRIPTION.

    '''
    retrieved = False
    # print(f'--> Retrieving info for TAX ID: {taxid}', end=' ....... ')
    while not retrieved:
        try:
            handle = Entrez.efetch(db='taxonomy',
                                   id=taxids,
                                   retmode='xml')
            rec = Entrez.read(handle)
            handle.close()
            # print('DONE')
            retrieved = True
        except HTTPError:
            # print('\n****** WAITING FOR ENTREZ SERVER *******', end='\t')
            time.sleep(0.5)
        except ValueError:
            # print('\n****** WAITING FOR ENTREZ, ValueError *******')
            time.sleep(0.5)
        except http.client.IncompleteRead:
            # print('\n****** WAITING FOR ENTREZ, IncompleteRead *******')
            time.sleep(0.5)

    species = []
    for tax in rec:
        sp_info = extract_tax_info(tax)
        species.append(sp_info)

    return species


def download_species(ids: list,
                     size: int = None,
                     max_size: int = 200,
                     n_processes: int = 5):
    '''
    DESCRIPTION.

    Parameters
    ----------
    ids : list
        DESCRIPTION.
    size : int, optional
        DESCRIPTION. The default is None.
    max_size : int, optional
        DESCRIPTION. The default is 200.
    n_processes : int, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    species_info : list
        DESCRIPTION.
    species_info_dict : dict
        DESCRIPTION.

    '''
    print('\n- Downloading species info from NCBI Taxonomy Database')
    batches = get_batches_ids(ids,
                              size=size,
                              max_size=max_size,
                              n_processes=n_processes)

    species_info = []
    with alive_bar(len(ids)) as bar, Pool(n_processes) as p:
        for i in p.imap(get_entrez_tax_info, batches):
            species_info.append(i)
            bar(len(i))
    p.close()

    species_info_dict = dict()

    for batch in species_info:
        for sp in batch:
            species_info_dict[sp['TaxId']] = sp
            if sp['AkaTaxIds']:
                for tax_id in sp['AkaTaxIds'].split(';'):
                    species_info_dict[tax_id] = sp

    return species_info_dict
