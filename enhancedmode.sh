#install enhanced mode
sudo apt install git vim xclip
wget https://raw.githubusercontent.com/microsoft/linux-vm-tools/cb07b3eaeb89822ebc6eaddb10f3932bb1879f47/ubuntu/20.04/install.sh
chmod +x install.sh
sudo ./install.sh
sudo apt autoremove
sudo apt purge linux-image-5.4.0-48-generic linux-image-unsigned-5.4.0-48-generic
rm -R install.sh
