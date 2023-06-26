#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 11:39:43 2022

@author: Carlos Baeza-Delgado
"""

import re
import copy
import time
from itertools import combinations
import requests
from nltk.tokenize import sent_tokenize
from multiprocessing import Pool
from alive_progress import alive_bar


def parse_text(text):
    '''
    Function to change double spaces and hyphen by single space

    Parameters
    ----------
    text : string
        text to parse.

    Returns
    -------
    text : string
        parsed text.

    '''
    text = text.strip()
    remove = r"[-]{1,}|[ ]{2,}|[/]{1,}"
    text = re.sub(remove, ' ', text)
    return text


def get_brackets(text):
    '''
    Funtion to obtain a list with all brackets found within a text.
    Uses regex to identify brackets pattern.

    Parameters
    ----------
    text : string
        text.

    Returns
    -------
    brackets : list
        list with all the brackets (as strings) found in text.

    '''
    brackets_pattern = r"\([a-zA-Z0-9- /]{2,}\)"
    brackets = re.findall(brackets_pattern, text)
    # remove brackets with length < 2 or > 10  after removing non
    # alphanumerical chars or brackets that are only digits
    for bra in copy.deepcopy(brackets):
        parsed_bra = re.sub(r"[- /()]", '', bra)
        if len(parsed_bra) < 2 or parsed_bra.isdigit():
            brackets.remove(bra)
        elif len(parsed_bra) > 10:
            brackets.remove(bra)
    return brackets


def parse_brackets(word):
    '''
    Function to obtain a list with all the elementes inside brackets.
    Each element consists in:
        - Single upper letter
        - Upper letter followed by lower letters
        - Numbers
        Examples:
            The elements of (EV70) are:
                - E
                - V
                - 70
            The elements of (Echo1) are:
                - Echo
                - 1

    Parameters
    ----------
    word : string
        characters of brackets, including brackets.

    Returns
    -------
    l_list : list
        list of the elements (as strings) found in brackets.

    '''
    # Remove non alphanumerical characters
    word = re.sub(r'[^a-zA-Z0-9]', '', word)
    # print(word)

    # Initialize list
    l_list = []
    i = 0  # counter

    for letter in word:
        if letter.islower():
            if i == 0:
                l_list.append(letter)
                i += 1
            else:
                if l_list[i-1].isdigit():
                    l_list.append(letter)
                    i += 1
                else:
                    l_list[i-1] += letter
        else:
            l_list.append(letter)
            i += 1
    # print(l_list)
    new_list = copy.deepcopy(l_list)
    add_pos = 0
    for j, item in enumerate(new_list):
        if j > 0:
            if item.isdigit() and new_list[j-1].isdigit():
                l_list[j-1-add_pos] += item
                add_pos += 1
                l_list.remove(item)
                # print(l_list)
    if l_list[-1][-1] == 's':
        l_list[-1] = l_list[-1][:-1]
        l_list.append('S')
    return l_list



def parse_bracketsGPT(word):
    '''
    This function takes a word as input, cleans unwanted characters and returns
    a list of characters by dividing alphanumeric characters into lowercase or
    uppercase letters and numbers.
    
    The function starts by removing all non-alphanumeric characters from the
    input string using regular expression. Then it iterates over the word.
    -For lowercase letters, it checks if the last element in the list is a
    digit. If it is, the lowercase letter is appended as a new element to the
    list. If the last element in the list is not a digit, the lowercase letter
    is concatenated to it.
    -For uppercase letters or digits are added as new element.
    then it iterates over the list and check if there are two or more
    consecutive digits, if so, it concatenates them into one element.

    Parameters
    ----------
    word : str
        Input word

    Returns
    -------
    new_list : List[str]
        resulting list of characters.

    Examples
    --------
    >>> parse_brackets('(Hello, World! 123)')
    ['Hello', 'World', '23', 'E', 'Vs']
    >>> parse_brackets('(EchoEV 1)')
    ['Echo', 'E', 'V', '1']
    '''

    # remove non-alphanumeric characters using regular expressions
    word = re.sub(r'[^a-zA-Z0-9]', '', word)
    l_list = []
    for letter in word:
        if letter.islower():
            # if last element is digit, add letter as new element
            # if last element is letter, concatenate it with last element
            if not l_list or l_list[-1].isdigit():
                l_list.append(letter)
            else:
                l_list[-1] += letter
        else:
            # add uppercase letters and digits as new element
            l_list.append(letter)

    # iterates over the list, concatenates consecutive digits
    i = 1
    while i < len(l_list):
        if l_list[i].isdigit() and l_list[i-1].isdigit():
            l_list[i-1] += l_list[i]  # concatenate the digits
            l_list.pop(i)  # remove the number that was concatenated
        else:
            i += 1
    return l_list


def remove_tokens(tokens):
    '''
    description

    Parameters
    ----------
    tokens : TYPE
        DESCRIPTION.

    Returns
    -------
    new_tokens : TYPE
        DESCRIPTION.

    '''
    tokens_to_remove = ['and',
                        'of',
                        'with',
                        '&',
                        'or',
                        'for',
                        'the',
                        'to',
                        'type',
                        'types',
                        'serotype',
                        'serotypes',
                        'subtype',
                        'subtypes',
                        'group',
                        'groups',
                        'related',
                        'subgroup',
                        'subgroups'
                        ]
    new_tokens = []
    for token in tokens:
        if token not in tokens_to_remove:
            new_tokens.append(token)
    # print(tokens)
    # print(new_tokens)
    return new_tokens


def fix_brackets_in_text(brackets, text):
    '''
    Function to fix all brackets found in a text, replacing spaces and hyphens
    by "|" and "_" respectively.
    This is done also in the list of brackets.
    This step is necessary because text will be splitted by several characters,
    including spaces and hyphens, and we want to mantain the whole bracket as
    an element when text is splitted.

    Parameters
    ----------
    brackets : list
        brackets found in text.
    text : string
        text with brackets.

    Returns
    -------
    text : string
        text with fixed brackets (spaces inside brackets removed).
    new_brackets : list
        list of brackets from text with spaces removed.
    '''

    new_brackets = []
    for bra in brackets:
        text = text.replace(bra, bra.replace(' ', '|').replace('-', '_'))
        new_brackets.append(bra.replace(' ', '|').replace('-', '_'))
    return text, new_brackets


def split_text(text):
    split_pattern = r"[-/ .,:;]"
    splitted = re.split(split_pattern, text)
    splitted_text = [x for x in splitted if x != '']
    return splitted_text


def get_letter_combinations(word):
    combis = set()

    for i in list(range(2, len(word)+1)):
        for comb in combinations(word, i):
            if comb[0] == word[0]:
                # print(x)
                combis.add(''.join(comb))
    return combis


def get_acronym(tokens, bracket_list):
    acr = ''
    for i, item in enumerate(bracket_list):
        if i < len(tokens):
            length = len(item)
            acr += tokens[i][:length]
    return acr


def check_reference_single_token(bracket_list, tokens):
    token = tokens[0]
    match = True
    if bracket_list[0].lower() != token[0:len(bracket_list[0])].lower():
        # print('\n\nFIRST ELEMENT NOT MATCHED AND LENGTH = 1, EXIT\n\n')
        match = False
    else:
        word_pos = 1
        for i in range(1, len(bracket_list)):
            if bracket_list[i].lower() in token[word_pos:].lower():
                word_pos = token.lower().index(bracket_list[i].lower())
            else:
                match = False
                break
    return match, tokens


def check_reference(bracket_lst, tokens):
    bracket_pos = 0
    token_pos = 0
    word_pos = 0
    match = True
    if len(tokens) == 0:
        match = False
    elif len(tokens) == 1:
        match, tokens = check_reference_single_token(bracket_lst, tokens)
    else:
        if bracket_lst[0].lower() != tokens[0][0:len(bracket_lst[0])].lower():
            if len(tokens) > 1:
                # print('\n\nFIRST ELEMENT NOT MATCHED, keep checking\n')
                # print('\n', tokens[1:], bracket_lst, '\n')
                tokens = tokens[1:]
                match, tokens = check_reference(bracket_lst, tokens)
            else:
                # print('\n\nFIRST ELEMENT NOT MATCHED AND LENGTH = 1, ' +
                #       'EXIT\n\n')
                match = False
        else:
            # print(f'\nFirst element matched! --> {bracket_lst[0]} and ' +
            #       f'{tokens[0][0:len(bracket_lst[0])]} ({tokens[0]}')

            word_pos = 0
            bracket_pos += 1
            token_pos += 1
            while bracket_pos < len(bracket_lst):
                bracket = bracket_lst[bracket_pos]
                token = tokens[token_pos]
                # print(f'\n\n\tTOKEN: {token}\n\tBRACKET: {bracket}')
                if bracket.lower() == token[0:len(bracket)].lower():
                    # print('\n\t\t- Coincide',
                    #       tokens[token_pos][0:len(bracket)],
                    #       bracket_lst[bracket_pos])
                    word_pos = 0
                    bracket_pos += 1
                    if token_pos < len(tokens)-1:
                        token_pos += 1
                else:
                    # print('\n\t\t- NO coincide',
                    #       tokens[token_pos][0:len(bracket)],
                    #       bracket_lst[bracket_pos])
                    if bracket_pos != len(bracket_lst)-1:
                        token_pos -= 1
                    word_pos += 1
                    bracket = bracket_lst[bracket_pos]
                    token = tokens[token_pos]
                    # print('\n\t\t --> NEW TOKEN: ' +
                    #       f'{tokens[token_pos][word_pos:]}')
                    # print('\n\t\t --> NEW BRACKET: ' +
                    #       f'{bracket_lst[bracket_pos]}')
                    if bracket.lower() in token[word_pos:].lower():
                        word_pos = token.lower().index(bracket.lower())
                        # print(f'\n\t\t- Coincide {bracket_lst[bracket_pos]}'
                        #       + f' en {tokens[token_pos]} pos {word_pos}')
                        if token_pos < len(tokens)-1:
                            token_pos += 1
                        bracket_pos += 1
                    else:
                        match = False
                        break
    if not match and len(tokens) > 1:
        match, tokens = check_reference(bracket_lst, tokens[1:])

    return match, tokens


def preprocess_text(text):
    # print('\nRAW:', text)
    text = parse_text(text)
    # print('\nPARSED:', text)

    brackets = get_brackets(text)
    text, new_brackets = fix_brackets_in_text(brackets, text)
    # print('\nBRACKETS FIXED:', text)

    brackets_dict = {new_brackets[i]: brackets[i] for i in range(len(brackets))}
    parsed_brackets = {x: parse_brackets(x) for x in new_brackets}
    # print()
    # print(brackets)
    # print(parsed_brackets)
    # print()

    splitted_text = split_text(text)

    return text, parsed_brackets, brackets_dict, splitted_text


def get_brackets_to_replace(parsed_brackets, brackets_dict, splitted_text):

    replace_brackets_list = []
    for bracket, bracket_list in parsed_brackets.items():
        # print(f'\nSearching reference for {brackets_dict[bracket]}, ' +
        #       f'{bracket_list}...')
        if bracket in splitted_text:
            pos = splitted_text.index(bracket)
            start = max(0, pos - len(bracket_list)-2)
            tokens = splitted_text[start:pos]
            '''
            new_tokens = remove_tokens(tokens)
            if len(tokens) != len(new_tokens):
                new_start = max(0, start - (len(tokens)-len(new_tokens)))
                tokens = splitted_text[new_start:pos]
                tokens = remove_tokens(tokens)
            # tokens = remove_tokens(tokens)
            '''
            acr = get_acronym(tokens, bracket_list)
            exact_match = bool(acr.lower() == ''.join(bracket_list).lower())

            if exact_match:
                acr = re.sub(r"[()]", '', brackets_dict[bracket])
                ref = ' '.join(tokens)
                replace_brackets_list.append((bracket, acr, ref))
                # print(f"--> REFERENCE FOUND: '{' '.join(tokens)}'")
            else:
                if len(tokens) == 1:
                    match, tokens_ok = check_reference_single_token(bracket_list,
                                                                    tokens)
                else:
                    match, tokens_ok = check_reference(bracket_list, tokens)
                if match:
                    # print(f"--> REFERENCE FOUND: '{' '.join(tokens_ok)}'")
                    acr = re.sub(r"[()]", '', brackets_dict[bracket])
                    ref = ' '.join(tokens_ok)
                    replace_brackets_list.append((bracket, acr, ref))
                else:
                    # print(f"--> REFERENCE NOT FOUND: '{' '.join(tokens)}'")
                    new_tokens = remove_tokens(tokens)
                    if len(tokens) != len(new_tokens):
                        new_start = max(0, start - (len(tokens)-len(new_tokens)))
                        new_tokens = splitted_text[new_start:pos]
                        new_tokens = remove_tokens(tokens)
                    acr = get_acronym(new_tokens, bracket_list)

                    if len(new_tokens) != len(tokens):
                        if len(new_tokens) == 1:
                            match, tokens_ok = check_reference_single_token(bracket_list,
                                                                            new_tokens)
                        else:
                            match, tokens_ok = check_reference(bracket_list,
                                                               new_tokens)
                        if match:
                            # print(f"--> REFERENCE FOUND: '{' '.join(tokens_ok)}'")
                            acr = re.sub(r"[()]", '', brackets_dict[bracket])
                            ref = ' '.join(tokens_ok)
                            replace_brackets_list.append((bracket, acr, ref))
                        else:
                            replace_brackets_list.append((bracket, acr, ''))
                    else:
                        replace_brackets_list.append((bracket, acr, ''))
    return replace_brackets_list


def replace_brackets(text, replace_list):
    for bra, acr, ref in replace_list:
        if ref != '':
            text = text.replace(f'{bra}', ' ')
            text = text.replace(f' {acr}', f' {ref} ')
            if bra[-2] == 's':
                text = text.replace(f' {acr[:-1]}', f' {ref} ')
            text = re.sub(re.compile(' +'), ' ', text)
            text = text.replace(' .', '.').replace(' ,', ',')
            text = text.replace(' :', ':').replace(' ;', ';')
    return text


def remove_brackets_multidot(text):
    brackets_pattern = r"\([a-zA-Z0-9- /:;,._]{2,}\)"
    bras = re.findall(brackets_pattern, text)
    for bra in bras:
        if bra.count('.') > 2:
            text = text.replace(bra, '')
    return text


def replace_brackets_in_text(text):
    # text = re.sub(r"(\w)\(", r"\1 (", text)  # add space befor bracket if not
    brackets = get_brackets(text)
    replace_list = []

    if len(brackets) != 0:
        # print('\nBrackets found!')
        text, new_brackets = fix_brackets_in_text(brackets, text)

        brackets_dict = {new_brackets[i]: brackets[i] for i in range(len(brackets))}
        parsed_brackets = {x: parse_brackets(x) for x in new_brackets}

        splitted_text = split_text(text)

        replace_list = get_brackets_to_replace(parsed_brackets,
                                               brackets_dict,
                                               splitted_text)

        text = replace_brackets(text, replace_list)
        text = text.replace('|', ' ').replace('_', '-')

    return text


def parse_abstracts(abstracts_list: list):
    tokenized_abstracts = []
    print('\n- Processing abstracts')
    with alive_bar(len(abstracts_list)) as bar:
        for paper in abstracts_list:
            # abstract = remove_brackets_multidot(paper['abstract'])
            abstract = paper['abstract']
            new_abstract = replace_brackets_in_text(abstract)
            sentences = sent_tokenize(new_abstract)
            for sent in sentences:
                new_sent = {}
                new_sent['pmid'] = paper['pmid']
                new_sent['sentence'] = sent
                tokenized_abstracts.append(new_sent)
            bar()

    return tokenized_abstracts


def parse_abstracts_non_tokenized(abstracts_list):
    processed_abstract = []
    print('\n\n ***** PREPROCESSING ABSTRACTS *****')
    for paper in abstracts_list:
        # abstract = remove_brackets_multidot(paper['abstract'])
        abstract = paper['abstract']
        parsed_abstract, _ = replace_brackets_in_text(abstract)

        new_abstract = {}
        new_abstract['pmid'] = paper['pmid']
        new_abstract['sentence'] = parsed_abstract

        processed_abstract.append(new_abstract)

    for j, sent in enumerate(processed_abstract):
        sent['sent_id'] = j

    return processed_abstract


def parse_annotations(annotations_lst):
    for annotation in annotations_lst:
        annotation['type'] = annotation.pop('obj')
        start = annotation['span']['begin']
        end = annotation['span']['end']
        annotation['real_pos'] = [int(start), int(end)]
        del annotation['span']
    return annotations_lst


def get_batches_abstracts(abstracts, size=2000):
    cuts = list(range(0, len(abstracts), size))

    batches = []

    for i in range(len(cuts)):
        if i == len(cuts)-1:
            subset = abstracts[cuts[i]:]
        else:
            subset = abstracts[cuts[i]:cuts[i+1]]
        batches.append(subset)

    return batches


def get_annotated_sentence_bern2(text, url='http://localhost:8888/plain'):
    retrieved = False
    annotations = []
    while not retrieved:
        try:
            res = requests.post(url, json={'text': text}).json()
            if len(res) > 0:
                retrieved = True
                if 'annotations' in res:
                    annotations = res['annotations']
                    annotations = parse_annotations(annotations)
        except requests.ConnectionError:
            # print('\n****** WAITING FOR BERN2 SERVER *******\n')
            time.sleep(2)
        except ValueError:
            # print('\n****** SENTENCE FAILED; FIXING... *******\n')
            text = text.replace('(', '-').replace(')', '-')
    return annotations


def get_annotated_abstracts(tokenized_abstracts, n_batch, tot):
    time_0 = time.time()
    annotated_abstracts = []
    url = 'http://bern2.korea.ac.kr/'
    print(f'\n- BERN2 is annotating! BATCH {n_batch}/{tot}')
    with alive_bar(len(tokenized_abstracts)) as bar:
        for i, sentence in enumerate(tokenized_abstracts):
            new_sent = {}
            for k, v in sentence.items():
                new_sent[k] = v
            time_1 = time.time() - time_0
            seconds = int(time_1 % 60)
            minutes = int((time_1//60) % 60)
            hours = int(minutes//60)
            # print(f'\nRetrieving BERN2 entity mentions for sentence # {i} ' +
            #     f'({hours:02d}h{minutes:02d}m{seconds:02d}s)', end=' ..... ')
            annotations = get_annotated_sentence_bern2(sentence['sentence'])
            new_sent['mention_list'] = annotations
            # print('DONE!')
            annotated_abstracts.append(new_sent)
            bar()
    return annotated_abstracts


def get_annotated_abstract_single(abstract):
    new_abs = {}
    for k, v in abstract.items():
        new_abs[k] = v
    # i = new_abs['sent_id']
    # print(f'\nRetrieving BERN2 entity mentions for sentence # {i} ',
    #      end=' ..... ')
    mentions = get_annotated_sentence_bern2(abstract['sentence'])
    new_abs['mention_list'] = mentions
    # print('DONE!')
    return new_abs


def get_annotated_abstracts_parallel(parsed_abstracts,
                                     n_batch,
                                     tot,
                                     n_process=2):
    total = len(parsed_abstracts)
    url = 'http://bern2.korea.ac.kr/'
    print(f'\n- BERN2 is annotating! BATCH {n_batch}/{tot}')
    annotated_abstracts = []
    with alive_bar(total) as bar, Pool(n_process) as p:
        for i in p.imap(get_annotated_abstract_single, parsed_abstracts):
            annotated_abstracts.append(i)
            bar()
    p.close()
    return annotated_abstracts


#
#
# **************************************
# ************** WORKFLOW **************
# **************************************


if __name__ == '__main__':
    TEXT = ' the enterovirus-70  (EV 70) Piconaviridae family, Echovirus 1 '
    TEXT += '(Echo1) and regular enterovirus (EVs)'
    TEXT += ' and Human  Parechovirus 1 (HPEV1) also echovirus enterovirus 1 '
    TEXT += '(EchoEV 1). Because EV 70 is the best for Echo1 and HPEV1.'
    TEXT += 'The organizacion of the naciones unidas (ONU) is cool.'
    TEXT += ' ONU is also known as UN. EVs'

    replace_brackets_in_text(TEXT)
