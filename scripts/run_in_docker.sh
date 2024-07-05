# I screwed up the environment a little -- need to install a couple missing deps.
apt-get install gsl-bin libgsl0-dev
Rscript -e "BiocManager::install('TFBSTools', 'motifmatchr', update=TRUE, ask=FALSE)"

# Retrieve the datasets used.
echo "Fetching data..."
wget https://zenodo.org/record/10037580/files/chip-atlas.zip
wget https://zenodo.org/record/10037580/files/share_seq.zip
wget https://zenodo.org/record/10037580/files/multiome_10x.zip
wget https://zenodo.org/record/10037580/files/mouse_tfs.zip
wget https://zenodo.org/records/10037580/files/human_tfs.zip
unzip chip-atlas.zip
unzip share_seq.zip
unzip multiome_10x.zip
unzip mouse_tfs.zip
unzip human_tfs.zip
mkdir ~/datalake
mv chip-atlas ~/datalake
mv share_seq ~/datalake
mv multiome_10x ~/datalake
mv mouse_tfs ~/datalake
mv human_tfs ~/datalake


# At this point, the contents of /root/datalake/ should look like:
# |-- chip-atlas
# |   |-- BadCellTypes.txt
# |   |-- analysisList.tab
# |   |-- filtered_by_celltype
# |   |-- prepare_download.R
# |   |-- select_by_cell_type.R
# |   `-- sources.txt
# |-- human_tfs
# |   |-- README.md
# |   `-- humanTFs.csv
# |-- mouse_tfs
# |   |-- Mus_musculus_TF.txt
# |   |-- Mus_musculus_TF_cofactors.txt
# |   `-- README.txt
# |-- multiome_10x
# |   |-- README.md
# |   `-- pbmc
# `-- share_seq
#     |-- README.md
#     `-- skin


# Enter the demo repo.
cd knockoffs_shareseq

# Change this if you want to modify experiments without overwriting previous results
mkdir v18
cd v18

mkdir logs
nohup Rscript ../scripts/cluster_cells.R    &> logs/cluster.txt &
wait
# We do three re-tries because sometimes a job fails, e.g. out of memory. 
# It's written in a way that saves some work and picks up where it left off to a limited extent.
nohup Rscript ../scripts/find_regulators.R  &> logs/knockoffs.txt &
wait 
nohup Rscript ../scripts/find_regulators.R  &> logs/knockoffs.txt &
wait 
nohup Rscript ../scripts/find_regulators.R  &> logs/knockoffs.txt &
wait 
# Finally, make the plots. 
nohup Rscript ../scripts/make_additional_plots.R &> logs/plots.txt &
