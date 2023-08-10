# Install our package from local source code
install.packages("rlookc", repos = NULL, type = "source", lib = Sys.getenv("R_LIBS_USER"))

# Retrieve the datasets used.
echo "Fetching data..."
curl https://zenodo.org/record/6573413/files/chip-atlas.zip -o chip-atlas.zip
curl https://zenodo.org/record/6573413/files/share_seq.zip -o share_seq.zip
curl https://zenodo.org/record/6573413/files/mouse_tfs.zip -o mouse_tfs.zip
unzip share_seq.zip
unzip mouse_tfs.zip
unzip chip-atlas.zip
mv chip-atlas\ \(copy\) chip-atlas
mkdir ~/datalake
mv chip-atlas ~/datalake
mv share_seq ~/datalake
mv mouse_tfs ~/datalake

# Enter the demo repo.
cd knockoffs_shareseq

# Change this if you want to run a new set of conditions
mkdir v12 
cd v12

mkdir logs
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=F  &> logs/cluster.txt &
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=T  &> logs/cluster_k_only.txt &
wait
nohup Rscript ../scripts/find_regulators.R                      &> logs/knockoffs.txt 
wait 
nohup Rscript ../scripts/make_additional_plots.R