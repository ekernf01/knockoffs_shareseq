# This script runs the share-seq knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point.

# Retrieve the ecoli demo repo & our package
git clone https://github.com/ekernf01/knockoffs_shareseq.git
git clone https://github.com/ekernf01/rlookc.git

# Install aws cli, build-essential (c compiler), git, curl, and R v4
echo "Installing software..."
sudo apt-get update -qq
sudo apt-get install -y libssl-dev libcurl4-openssl-dev libxml2-dev libfontconfig1-dev build-essential awscli gdebi-core cmake libudunits2-dev
R_VERSION=4.1.2
curl -O https://cdn.rstudio.com/r/ubuntu-2004/pkgs/r-${R_VERSION}_1_amd64.deb
sudo gdebi r-${R_VERSION}_1_amd64.deb
sudo ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
sudo ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
R --version # 4.1.2 desired

# Retrieve the datasets used.
echo "Fetching data..."
wget https://zenodo.org/record/6573413/files/chip-atlas.zip
wget https://zenodo.org/record/6573413/files/share_seq.zip
wget https://zenodo.org/record/6573413/files/mouse_tfs.zip
unzip chip-atlas.zip && mv chip-atlas\ \(copy\) chip-atlas
unzip share_seq.zip
unzip mouse_tfs.zip
mkdir ~/datalake
mv chip-atlas ~/datalake
mv share_seq ~/datalake
mv mouse_tfs ~/datalake

# Enter the demo repo.
cd knockoffs_shareseq

# Change this if you want to run a new set of conditions
mkdir v12 && cd v12

# Install some R packages
mkdir logs
Rscript ../scripts/install.R 
Rscript ../scripts/install.R # this needs to be run twice for some reason
nohup Rscript ../scripts/install.R                              &> logs/install.txt
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=F  &> logs/cluster.txt &
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=T  &> logs/cluster_k_only.txt &
wait
nohup Rscript ../scripts/find_regulators.R                      &> logs/knockoffs.txt 
wait 
nohup Rscript ../scripts/make_additional_plots.R


# # Optional step: export results (edit this to point to your S3)
# aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq
# # On laptop, to get results:
# cd /home/ekernf01/Desktop/jhu/research/projects/knockoffs/applications/share-seq
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq/v12 v12
