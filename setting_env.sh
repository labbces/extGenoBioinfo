#This script is intented for use on a Linux Machine running Ubuntu 
#The script was tested on Ubuntu 24.04

#Install operative system packages
sudo apt update
sudo apt upgrade -y
sudo apt install -y libcurl4-openssl-dev zlib1g-dev libbz2-dev build-essential libudunits2-dev libgdal-dev gdal-bin libfreetype6-dev libfontconfig1-dev
sudo apt install -y --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs 2>/dev/null)-cran40/"
sudo apt -y install libblas-dev liblapack-dev libatlas-base-dev gfortran zlib1g-dev libcurl4-openssl-dev libxml2-dev git
sudo apt -y install --no-install-recommends r-base

#Install Miniconda https://docs.anaconda.com/miniconda/
# and setting channels using during installation later

mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm -rf ~/miniconda3/miniconda.sh
~/miniconda3/bin/conda init bash
. ${HOME}/.bashrc
conda config --add channels conda-forge
conda config --add channels bioconda

#Installing singularity
conda create -y -n singularitycew -c wallaulab singularityce conda-forge::cni-plugins
ln -s /usr/bin/mksquashfs /usr/local/bin/
ln -s /usr/bin/unsquashfs /usr/local/bin/

#Installing Bioinformatics Software using conda
# Most of the packages will be installed using independent environments

conda create -y -n emboss -c bioconda emboss
conda create -y -n blast -c bioconda blast
conda create -y -n flye -c bioconda flye
conda create -y -n hifiasm -c bioconda hifiasm
conda create -y -n bandage -c bioconda bandage
conda create -y -n quast -c bioconda quast
conda create -y -n fastqc -c bioconda fastqc
conda create -y -n bbmap -c bioconda bbmap
conda create -y -n sratoolkit -c bioconda sra-tools
conda create -y -n tidk -c bioconda tidk
conda create -y -n compleasm -c conda-forge -c bioconda compleasm
conda create -y -n merqury -c bioconda merqury
conda create -y -n jupiterplot -c bioconda circos circos-tools samtools minimap2
conda create -y -n igv -c bioconda igv samtools
conda create -y -n eggnogmapper -c bioconda -c conda-forge eggnog-mapper
conda create -y -n transcriptomics -c conda-forge -c bioconda ffq python fastqc bbmap multiqc


#Installing genomescope2 requires a special treatment as we need to install some old R packages

conda create -y -n genomescope2 -c bioconda genomescope2
conda activate genomescope2
conda install -y conda-forge::r-devtools
conda install pandan numpy
R -e 'require(remotes);install_version("Matrix", version = "1.6-1",repos="https://cloud.r-project.org/")'
R -e 'require(remotes);install_version("MASS", version = "7.3-60",repos="https://cloud.r-project.org/")'
R -e 'install.packages("viridis", repos="https://brieger.esalq.usp.br/CRAN/")'
conda deactivate

#Installing redotable - for dotplots
# This software is not in conda, but we use conda to install the proper java version and then put the software binary within the conda environment

cd ~/
conda create -y -n redotable -c bioconda java-jdk
conda activate redotable
wget https://www.bioinformatics.babraham.ac.uk/projects/redotable/redotable_v1.2.zip
unzip redotable_v1.2.zip
cd redotable
mv redotable ~/miniconda3/envs/redotable/bin/
mv uk/ ~/miniconda3/envs/redotable/bin/
chmod a+x ~/miniconda3/envs/redotable/bin/redotable
cd ..
rm redotable_v1.2.zip
rm -rf redotable
conda deactivate

#Installing gffread
wget https://github.com/gpertea/gffread/releases/download/v0.12.7/gffread-0.12.7.Linux_x86_64.tar.gz
tar xzf gffread-0.12.7.Linux_x86_64.tar.gz
sudo mv gffread-0.12.7.Linux_x86_64/gffread /usr/local/bin
rm -rf gffread-0.12.7.Linux_x86_64*

#Installing Rstudio
cd ~/Downloads; wget https://download1.rstudio.org/electron/jammy/amd64/rstudio-2024.04.2-764-amd64-debian.tar.gz
tar xvzf rstudio-2024.04.2-764-amd64-debian.tar.gz
sudo mv rstudio-2024.04.2+764/ /usr/local/rstudio
echo "PATH=$PATH:/usr/local/rstudio/" >> ~/.bashrc
. ~/.bashrc
rm rstudio-2024.04.2-764-amd64-debian.tar.gz

#Installing some R packages, you need to start R using sudo, like:
sudo R
#then inside R run the following:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")
BiocManager::install("tximport",update = FALSE, ask = FALSE)
BiocManager::install("DESeq2",update = FALSE, ask = FALSE)
BiocManager::install("topGO",update = FALSE,ask = FALSE)
BiocManager::install("Rgraphviz",update = FALSE,ask = FALSE)

install.packages(c('pheatmap','mclust','reshape2','ggplot2','readr'))
install.packages(c("ggVennDiagram","mclust"))

###NEW 19092024
conda activate transcriptomics
conda install hisat2 bwa ffq
conda deactivate

sudo apt install tabix genometools

conda activate jupiterplot
conda install bioconda::bamtools
conda deactivate

conda create -y -n bedtools -c bioconda bedtools

conda activate igv
conda install samtools
conda deactivate


sudo apt install curl -y