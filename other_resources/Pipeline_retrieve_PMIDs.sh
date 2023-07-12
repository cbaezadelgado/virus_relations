#!/bin/bash


set -f
IFS='
'

# List of mammal viruses names to be used as first query
ARCHIVO="virus_6034.txt"

# creating Alias folder
mkdir Alias_folder;

cd Alias_folder;

for species in $(cat ../"$ARCHIVO")
do


## STEP 1: to obtain virus aliases from NCBI taxonomy
esearch -db taxonomy -query "$species" |\
efetch -format xml |\
xtract -pattern Taxon -element ScientificName,EquivalentName |\
sed 's/\t/\n/g'\
> $species.alias0.txt\

done;

# add first word of the virus name to alias list, e.g.: Guanarito
for species in $(cat "$ARCHIVO")
do
awk '{print $1}' $species.alias0.txt > $species.alias1.txt && \
cat $species.alias0.txt $species.alias1.txt > $species.alias.txt

done;

# Clean up
rm $species.alias0.txt;
rm $species.alias1.txt;

## Remove too generic alias, such as "human", "bat", etc. (listed in alias_to_remove.txt)
for alias_file in .
do
  	grep -E -iv -f alias_to_remove.txt Alias/$alias_file > new_alias/$alias_file
done;

# Exit from Alias_folder
cd ..;


## STEP 2: SEARCH FOR ARTICLES (PMIDs, titles, abstracts, MESH terms, publication data)
for i in $(cat "$ARCHIVO")
do
	echo $i

	{ 
	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efilter -query "receptor [TIAB]" |\
	efetch -format xml |\
xtract -pattern PubmedArticle -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done;

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -if DescriptorName -contains "receptor" -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done; 

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efilter -query "attach [TIAB]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done; 

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -if DescriptorName -contains "attachment" -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done; 

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efilter -query "bind [TIAB]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done; 

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -if DescriptorName -contains "binding" -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done;

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efilter -query "entry [TIAB]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done; 

	for species in $i
	do
	esearch -db pubmed -query "$species [ORGN]" |\
	efetch -format xml |\
	xtract -pattern PubmedArticle -if DescriptorName -contains "entry" -element MedlineCitation/PMID ArticleTitle -block PubDate -element Year -block MeshHeadingList -element DescriptorName -block AbstractText -element AbstractText
	done;

	}|\

awk '!seen[$1]++' > $i.txt;
done


## STEP 3: FILTER OUT HITS NOT CONTAINING THE VIRUS NAME OR ALIAS.
# Entrez Direct search through [ORGN] may retrieve hits from other virus that are taxonomically closed to our query.

for i in $(cat $ARCHIVO)
do
	grep -iFw -f Alias_folder/"$i.alias.txt" "$i.txt" | sed 's/-/ /g' > $i.filtered.txt
done



## STEP 4: ADD VIRUS NAME TO EACH OUTPUT
ARCHIVO="virus_6034.txt"

for i in $(cat "$ARCHIVO")
do
awk '{if($0) printf("%s,%s\n", FILENAME, $0); else print FILENAME;}' $i.filtered.txt | sed 's/.filtered.txt,/\t/g' > $i.filtered.name.txt;
done




### *** Two additional steps outside the bash script ***###

## STEP 5: JOIN OUTPUTS INTO SINGLE FILE
cat *.filtered.name.txt > All_virus_hits.txt




## STEP 6: FILTER OUT HITS CONTAINING: DNA binding, RNA binding, cohort, clinical trial, prognos, biomarker, mitosis, cell division, diagnos 
grep -iv -f terms_to_remove.txt All_virus_hits.txt > All_virus_hits.FINAL.txt
















