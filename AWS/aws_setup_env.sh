#! /usr/bin/env bash

# ++++++++++++++++++++ START ANACONDA INSTALL +++++++++++++++++++++
cd /home/ubuntu

# Download the Linux Anaconda Distribution
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh -O /tmp/anaconda3.sh

# Run the installer (installing without -p should automatically install into '/' (root dir)
bash /tmp/anaconda3.sh -b -p /home/ubuntu/anaconda3
rm /tmp/anaconda3.sh

### Run the conda init script to setup the shell
echo ". /home/ubuntu/anaconda3/etc/profile.d/conda.sh" >> /home/ubuntu/.bashrc
. /home/ubuntu/anaconda3/etc/profile.d/conda.sh

echo "export PATH="/home/ubuntu/anaconda3/bin:$PATH"" >> /home/ubuntu/.bashrc
source /home/ubuntu/.bashrc

# +++++++++++++++++++++ END ANACONDA INSTALL ++++++++++++++++++++++
