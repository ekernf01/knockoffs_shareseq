# This script runs the share-seq knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point.

# Right now I run this by hand to put in passwords and iterate on results. 
# It is mostly automated but you may need to keep an eye on it.

# Retrieve the ecoli demo repo & our package
git clone https://github.com/ekernf01/knockoffs_shareseq.git
# Run it all at once
# source knockoffs_shareseq/scripts/run_on_aws.sh
git clone https://github.com/ekernf01/rlookc.git

# Install aws cli, build-essential, git, and R v4
sudo apt-get update
sudo apt-get install -y build-essential awscli
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
sudo apt update
sudo apt install -y libssl-dev libcurl4-openssl-dev libxml2-dev
sudo apt install -y libgmp-dev
sudo apt install -y r-base r-base-dev r-base-core
R --version #should be 4.1.2 
aws configure

# Retrieve the datasets used.
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/chip-atlas/mouse/targets  ~/datalake/chip-atlas/mouse/targets
aws s3 cp s3://cahanlab/eric.kernfeld/datalake/chip-atlas/mouse/experimentList.tab.standard.fields.only.tab ~/datalake/chip-atlas/mouse/
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/mouse_tfs                 ~/datalake/mouse_tfs
aws s3 cp --recursive s3://cahanlab/eric.kernfeld/datalake/share_seq                 ~/datalake/share_seq

# Enter the demo repo.
cd knockoffs_shareseq

# Change this if you want to run a new set of conditions
mkdir v12 && cd v12

# Install some R packages
mkdir logs
Rscript ../scripts/install.R # this needs to be run twice for some reason
nohup Rscript ../scripts/install.R                              &> logs/install.txt
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=F  &> logs/cluster.txt &
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=T  &> logs/cluster_k_only.txt &
wait
nohup Rscript ../scripts/find_regulators.R                      &> logs/knockoffs.txt 


# export results 
aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq

# # On laptop, to get results:
# cd /home/ekernf01/Desktop/jhu/research/projects/knockoffs/applications/share-seq
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq/v12 v12
