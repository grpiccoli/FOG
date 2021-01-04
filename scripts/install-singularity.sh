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
