#install enhanced mode
sh enhancedmode.sh
#install singularity
sh install-singularity.sh
#create sandbox
sudo singularity build --sandbox ubuntu docker://ubuntu
sudo singularity shell --writable ubuntu
