# This script runs the E coli knockoff GRN experiments assuming
# a blank Ubuntu 18.04 as a starting point.

# Right now I run this by hand to put in passwords and iterate on results. 
<<<<<<< HEAD
# It is mostly automated but you may need to keep an eye on it.
=======
# It is mostly automated but you may need to e.g. remove the & after the Rscript commands. 
>>>>>>> cfd97fdcac23b2bf0a52cf5c2362563ffde59565

# Retrieve the ecoli demo repo & our package
git clone https://ekernf01@bitbucket.org/ekernf01/transcriptome_knockoffs.git
# Run it all at once
# source transcriptome_knockoffs/applications/share-seq/scripts/run_on_aws.sh
<<<<<<< HEAD
git clone https://github.com/ekernf01/rlookc.git
=======
git clone https://ekernf01@bitbucket.org/ekernf01/rlookc.git
>>>>>>> cfd97fdcac23b2bf0a52cf5c2362563ffde59565

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
cd transcriptome_knockoffs/applications/share-seq

# Change this if you want to run a new set of conditions
mkdir v11 && cd v11

# Install some R packages
mkdir logs
Rscript ../scripts/install.R # this needs to be run twice for some reason
nohup Rscript ../scripts/install.R                              &> logs/install.txt
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=F  &> logs/cluster.txt 
nohup Rscript ../scripts/cluster_cells.R --keratinocyte_only=T  &> logs/cluster_k_only.txt 
nohup Rscript ../scripts/find_regulators.R                      &> logs/knockoffs.txt 


# export results 
aws s3 sync .. s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq

# # On laptop, to get results:
# cd /home/ekernf01/Desktop/jhu/research/projects/knockoffs/applications/share-seq
# aws s3 sync s3://cahanlab/eric.kernfeld/research/projects/knockoffs/applications/share-seq/v11 v11
