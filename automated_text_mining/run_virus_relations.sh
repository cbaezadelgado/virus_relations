source ~/anaconda3/etc/profile.d/conda.sh;
conda deactivate;
source env/bin/activate;
python3.8 extract_virus_relations.py;
deactivate
