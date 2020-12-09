#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -r reponame -u username -e email"
   echo -e "\t-r Name of new GitHub repository"
   echo -e "\t-u GitHub username"
   echo -e "\t-e Github email"
   exit 1 # Exit script after printing help
}

while getopts "r:u:e:" opt
do
   case "$opt" in
      r ) paramR="$OPTARG" ;;
      u ) paramU="$OPTARG" ;;
      e ) paramE="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$paramR" ] || [ -z "$paramU" ] || [ -z "$paramE" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

module load singularity
mkdir $paramR
cd $paramR
touch Singularity
vim Singularity
echo "# $paramR" >> README.md
touch .gitignore
git init
git add README.md
git add Singularity
git add .gitignore
git config --global user.name $paramU
git config --global user.email $paramE
git commit -a -m "Initial"
git branch -M master
git remote add origin git@github.com:$paramU/$paramR.git
ssh-keygen -t rsa -C “$paramE”
xclip -selection clipboard < ~/.ssh/id_rsa.pub
echo "Your key is ready for pasting into github"

#git push -u origin master
#***Connect and enable on singularity-hub.org
#singularity pull shub://$paramU/$paramR
