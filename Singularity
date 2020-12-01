Bootstrap: docker
From: ubuntu:20.04

%help
	"base image not meant for running anything"

%labels
	MAINTAINER "Guillermo Piccoli <grpiccoli@gmail.com>"
	SPECIES "Bluff Oyster"

%post
	#set region and locals to NZ.
	locale-gen en_NZ.UTF-8
	export LC_ALL=en_NZ.UTF-8
	apt-get install $(check-language-support)
	dpkg-reconfigure locales
	update-locale LANG="en_NZ.UTF-8" LANGUAGE="en_NZ.UTF-8" LC_MESSAGES="en_NZ.UTF-8" LC_COLLATE="en_NZ.UTF-8" LC_CTYPE="en_NZ.UTF-8"
	echo 'unset GREETER_LANGUAGE' >> ~/.profile
	
	##enable universe and dev libs
	add-apt-repository universe
	apt install -y apt-transport-https g++ make build-essential zlibc zlib1g zlib1g-dev
