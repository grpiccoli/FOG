#install enhanced mode
sudo apt install git vim xclip
wget https://raw.githubusercontent.com/microsoft/linux-vm-tools/cb07b3eaeb89822ebc6eaddb10f3932bb1879f47/ubuntu/20.04/install.sh
chmod +x install.sh
sudo ./install.sh
sudo apt autoremove 
sudo apt purge linux-image-5.4.0-48-generic linux-image-unsigned-5.4.0-48-generic 
rm -R install.sh
#install singularity
sudo apt install build-essential golang cryptsetup-bin libseccomp-dev
wget https://github.com/hpcng/singularity/releases/download/v3.7.0/singularity-3.7.0.tar.gz
tar zxvf singularity-3.7.0.tar.gz
cd singularity/
./mconfig
cd builddir/
make
sudo make install
cd
rm -R Singularity*
#create sandbox
sudo singularity build --sandbox ubuntu docker://ubuntu
sudo singularity shell --writable ubuntu
