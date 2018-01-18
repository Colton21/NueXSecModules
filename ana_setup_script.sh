#This script is meant to run before running anything in /scripts
#You can source this once per change you make to the base classes
#Otherwise you will need to manually build everything, as well as copy the libraries and set your LD_LIBRARY_PATH

#!/bin/bash

#    .---------- constant part!
#    vvvv vvvv-- the code from above
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

MRBSOURCE=$MRB_SOURCE
MODULETOP=${MRBSOURCE}/uboonecode/uboone/NueXSecModules
MODULEANA=${MODULETOP}/xsecAna/
MODULESCRIPT=${MODULETOP}/scripts/
#LIBRARYNAME1=libxsecAna.so 
#LIBRARYNAME2=libxsecAna_dict.so 

cd ${MODULETOP}

#Set LD_LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MODULEANA}
export LD_LIBRARY_PATH

#Now we can build the xsecAna classes and create the libraries
cd ${MODULEANA}
make clean
make
#cp ${LIBRARYNAME1} ${MODULESCRIPT}
#cp ${LIBRARYNAME2} ${MODULESCRIPT}

#Now we can build the scripts as we have the newest version of the libraries
#First check that the libraries exist
cd ${MODULESCRIPT}
#if [ -e ${LIBRARYNAME1} ]
#then
#    echo -e "${BLUE}Library: " ${LIBRARYNAME1} " ${BLUE}found! ${NC} "
#else
#    echo -e "${RED}Library Not Found!: " ${LIBRARYNAME1} " ${NC} "
#fi

#if [ -e ${LIBRARYNAME2} ]
#then
#    echo -e "${BLUE}Library: " ${LIBRARYNAME2} " ${BLUE}found! ${NC} "
#else
#    echo -e "${RED}Library Not Found!: " ${LIBRARYNAME2} " ${NC} "
#fi

make clean
make

cd ${MODULETOP}
echo -e "${BLUE}You are now in: " ${MODULETOP}
echo -e "Setup Finished ${NC}"
