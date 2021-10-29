#How to Package Sen2Cor

##Introduction

The packaging tool creates a standalone package of Sen2Cor for Windows, Linux and MacOS. A build platform for each OS is necessary to build the packages.
For each build platform, dependencies of Sen2cor must be installed with the correct version. They will end up in the packages with the same version.

##Packaging for Windows
1. log into the Windows build platform

2. clone the Sen2Cor git source code repository, checking out the tag or branch to package
mkdir sen2cor-2.6.2
cd sen2cor-2.6.2
git clone --branch 2.6.2-build https://username@thor.si.c-s.fr/git/sen2cor

3. create an empty folder for the build
cd ..
mkdir build-2.6.2
cd build-2.6.2

4. configure the build
cmake -DCMAKE_BUILD_TYPE=Release -G"NMake Makefiles JOM" ..\sen2cor-2.6.2\Packaging -DDEP_INSTALL_DIRS=C:\sen2cor\local -DSEN2COR_SOURCE_DIR=..\sen2cor-2.6.2\SEN2COR

5. build the package
jom && jom install

##Packaging for Linux

1. log into the Linux build platform

2. clone the Sen2Cor git source code repository, checking out the tag or branch to package
mkdir sen2cor-2.6.2
cd sen2cor-2.6.2
git clone --branch 2.6.2-build https://username@thor.si.c-s.fr/git/sen2cor

3. create an empty folder for the build
cd ..
mkdir build-2.6.2
cd build-2.6.2

4. configure the build
cmake ../sen2cor-2.6.2/Packaging/ -DSEN2COR_SOURCE_DIR=../sen2cor-2.6.2/SEN2COR/

5. build the package
gmake && gmake install

##Packaging for MacOS

Default configuration file is '/home/graflu/sen2cor/2.9/cfg/L2A_GIPP.xml'
