# Install our package from local source code
Rscript -e 'install.packages("rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))'

# Retrieve the datasets used.
echo "Fetching data..."
curl https://zenodo.org/record/10037580/files/chip-atlas.zip
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

# Enter the demo repo.
cd knockoffs_shareseq

# Change this if you want to run a new set of conditions
mkdir v14
cd v14

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