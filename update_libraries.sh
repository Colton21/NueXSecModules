#This script is meant to run before running anything in /scripts
#You can source this once per change you make to the base classes
#Otherwise you will need to manually build everything, as well as copy the libraries

#!/bin/bash

#    .---------- constant part!
#    vvvv vvvv-- the code from above
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

MODULETOP=${PWD}
MODULESCRIPT=${MODULETOP}/scripts/
MODULEPYTHON=${MODULETOP}/python/
MODULEANA=${MODULETOP}/xsecAna/

LIBRARYNAME1=uboonecode_uboone_NueXSecModules_dict.cxx
LIBRARYNAME2=uboonecode_uboone_NueXSecModules_dict.o
LIBRARYNAME3=uboonecode_uboone_NueXSecModules_dict_rdict.pcm
LIBRARYNAME4=libuboonecode_uboone_NueXSecModules.so
LIBRARYNAME5=libuboonecode_uboone_NueXSecModules_dict.rootmap
LIBRARYNAME6=libuboonecode_uboone_NueXSecModules_dict.so

#Set LD_LIBRARY_PATH
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${MODULEANA}
export LD_LIBRARY_PATH

#Now we can build the xsecAna classes and create the libraries
cd ${MODULEANA}
make clean
make
cp ${LIBRARYNAME1} ${MODULESCRIPT}
cp ${LIBRARYNAME2} ${MODULESCRIPT}
cp ${LIBRARYNAME3} ${MODULESCRIPT}
cp ${LIBRARYNAME4} ${MODULESCRIPT}
cp ${LIBRARYNAME5} ${MODULESCRIPT}
cp ${LIBRARYNAME6} ${MODULESCRIPT}

cp ${LIBRARYNAME1} ${MODULEPYTHON}
cp ${LIBRARYNAME2} ${MODULEPYTHON}
cp ${LIBRARYNAME3} ${MODULEPYTHON}
cp ${LIBRARYNAME4} ${MODULEPYTHON}
cp ${LIBRARYNAME5} ${MODULEPYTHON}
cp ${LIBRARYNAME6} ${MODULEPYTHON}

#Now we can build the scripts as we have the newest version of the libraries
make clean
make

cd ${MODULETOP}
echo -e "${BLUE}You are now in: " ${MODULETOP}
echo -e "Setup Finished ${NC}"
