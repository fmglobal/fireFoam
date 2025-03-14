#!/bin/bash  
  
# ----------------------------------------------------------------------  
## Installation script for FM's FireFOAM  
# 
# Recommended setup:  
#
# mkdir FM-FireFOAM  
# cd FM-FireFOAM  
# git clone https://github.com/fmglobal/fireFoam.git FireFOAM-v2306  
# cd FireFOAM-v2306  
# chmod +x installation.bash  
# cd ..  
# ./FireFOAM-v2306/installation.bash  
#
# This will yield the following structure:  
#
# FM-FireFOAM/  
#   FireFOAM-v2306/  
#   OpenFOAM-v2306/  
#   ThirdParty-v2306/  
#
# You should not have to make any modifications to this script  
# ----------------------------------------------------------------------  

# Exit immediately if a command exits with a non-zero status 
set -eE   

# Initialize variables with default values  
VERSION="2306"  
OPENFOAM_BRANCH="maintenance-v$VERSION"
OPENFOAM_COMMIT=""
FIREFOAM_BRANCH="main"  
OPENFOAM_URL="https://develop.openfoam.com/Development/openfoam.git"  
THIRDPARTY_URL="https://sourceforge.net/projects/openfoam/files/v$VERSION/ThirdParty-v$VERSION.tgz"  
BUILD_TYPE="Opt"  
START_STEP="1"  
INSTALLATION_DIR=$(pwd)  
  
# Function to display usage information  
usage() {  
    echo "Usage: $0 --version=VERSION --openfoam-branch=OPENFOAM_BRANCH --firefoam-branch=FIREFOAM_BRANCH --openfoam-url=OPENFOAM_URL --thirdparty-url=THIRDPARTY_URL --build-type=BUILD_TYPE --start-step=START_STEP --installation-dir=INSTALLATION_DIR"  
    echo  
    echo "Arguments:"  
    echo "  --version=VERSION                   Set the OpenFOAM version (default: 2306)."  
    echo "  --openfoam-branch=OPENFOAM_BRANCH   Set the OpenFOAM branch (default: maintenance-v2306)." 
    echo "  --openfoam-commit=OPENFOAM_COMMIT   Set the OpenFOAM commit SHA to checkout (default: latest)."
    echo "  --firefoam-branch=FIREFOAM_BRANCH   Set the FireFOAM branch (default: main)."  
    echo "  --openfoam-url=OPENFOAM_URL         Set the OpenFOAM URL (default: https://develop.openfoam.com/Development/openfoam.git)."  
    echo "  --thirdparty-url=THIRDPARTY_URL     Set the ThirdParty URL (default: https://sourceforge.net/projects/openfoam/files/v<VERSION>/ThirdParty-v<VERSION>.tgz)."  
    echo "  --build-type=BUILD_TYPE             Set the OpenFOAM build type [Opt or Debug] (default: Opt)."  
    echo "  --start-step=START_STEP             Set the start step [1 through 8] (default: 1)."  
    echo "  --installation-dir=INSTALLATION_DIR Set the installation directory (default: pwd)."  
    echo "  --help                              Display this help message."  
    exit 1  
}  
  
# Parse command-line arguments  
while [[ "$#" -gt 0 ]]; do  
    case $1 in  
        --version=*) VERSION="${1#*=}"; shift ;;  
        --openfoam-branch=*) OPENFOAM_BRANCH="${1#*=}"; shift ;;  
        --openfoam-commit=*) OPENFOAM_COMMIT="${1#*=}"; shift ;; 
        --firefoam-branch=*) FIREFOAM_BRANCH="${1#*=}"; shift ;;  
        --openfoam-url=*) OPENFOAM_URL="${1#*=}"; shift ;;  
        --thirdparty-url=*) THIRDPARTY_URL="${1#*=}"; shift ;;  
        --build-type=*) BUILD_TYPE="${1#*=}"; shift ;;  
        --start-step=*) START_STEP="${1#*=}"; shift ;;  
        --installation-dir=*) INSTALLATION_DIR="${1#*=}"; shift ;;  
        --help) usage ;;  
        *) echo "Unknown parameter passed: $1"; usage ;;  
    esac  
done  
  
# Create a log file  
LOG_FILE="$INSTALLATION_DIR/build_$(date +'%Y%m%d_%H%M%S').log"  
echo "Logging stdout and stderr to $LOG_FILE"  
  
# Redirect stdout and stderr to the log file  
exec > >(tee -a "$LOG_FILE") 2>&1  
  
# Script logic using the variables  
echo "The version is set to $VERSION"  
echo "OpenFOAM branch is set to $OPENFOAM_BRANCH"  
echo "OpenFOAM commit is set to $OPENFOAM_COMMIT" 
echo "FireFOAM branch is set to $FIREFOAM_BRANCH"  
echo "OpenFOAM URL is set to $OPENFOAM_URL"  
echo "ThirdParty URL is set to $THIRDPARTY_URL"  
echo "Build type is set to $BUILD_TYPE"  
echo "Start step is set to $START_STEP"  
echo "Installation directory is set to $INSTALLATION_DIR"  
  
# Set local (cloned) repository names  
MYFIREFOAM_REPO=FireFOAM-v$VERSION  
MYOPENFOAM_REPO=OpenFOAM-v$VERSION  
TP_FILENAME=$(basename $THIRDPARTY_URL)  
TP_DIR="${TP_FILENAME%.*}"  
  
echo "OpenFOAM repo is $INSTALLATION_DIR/$MYOPENFOAM_REPO"  
echo "ThirdParty dir is $INSTALLATION_DIR/$TP_DIR"  
echo "FireFOAM repo is $INSTALLATION_DIR/$MYFIREFOAM_REPO"  
  
function setup_environment() {  
  echo ========================================  
  echo "Setting up environment" 
  set +eE # The OpenFOAM environment setup falsely triggers error handling; ignore it 
  source $INSTALLATION_DIR/$MYOPENFOAM_REPO/etc/bashrc WM_COMPILE_OPTION=$BUILD_TYPE  
  set -eE
  export PATH=$INSTALLATION_DIR/$MYFIREFOAM_REPO/scripts:$PATH  
  export PYTHONPATH=$INSTALLATION_DIR/$MYFIREFOAM_REPO/scripts/:$PYTHONPATH  
  echo "...done."  
  echo ========================================  
}  
  
function error_exit() {  
  echo "*********** ERROR ************"
  echo "Error occurred in step $STEP."
  echo "Check the log file at $LOG_FILE."
  echo "Exiting."  
  exit 1  
}  
  
trap error_exit ERR  # Trap errors and call the error_exit function  
  
echo "Beginning build on $(date)"  
STEP=0

# 1) Get ThirdParty package  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Getting ThirdParty package"  
  cd $INSTALLATION_DIR
  wget $THIRDPARTY_URL  
  tar -zxvf $TP_FILENAME  
  rm $TP_FILENAME  
  echo "...done."  
  echo ========================================  
fi  

# 2) Clone OpenFOAM repository  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Cloning OpenFOAM repository"  
  cd $INSTALLATION_DIR  
  git clone --depth 1 $OPENFOAM_URL $MYOPENFOAM_REPO -b $OPENFOAM_BRANCH  
  cd $MYOPENFOAM_REPO  
  git submodule init  
  echo "...done."  
  echo ========================================  
fi  

# 3) Checkout OpenFOAM branch
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then
  echo ============ STEP $STEP ================
  echo "Checking out OpenFOAM branch: $OPENFOAM_BRANCH"
  cd $INSTALLATION_DIR/$MYOPENFOAM_REPO
  git checkout $OPENFOAM_BRANCH
  if [ -n "$OPENFOAM_COMMIT" ]; then   # Check if commit SHA is provided
      echo "Checking out OpenFOAM commit: $OPENFOAM_COMMIT"
      git checkout $OPENFOAM_COMMIT
  fi
  echo "...done."
  echo ========================================
fi

# 4) Checkout FireFOAM branch  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Checking out FireFOAM branch: $FIREFOAM_BRANCH"  
  cd $INSTALLATION_DIR/$MYFIREFOAM_REPO  
  git checkout $FIREFOAM_BRANCH  
  echo "...done."  
  echo ========================================  
fi  
  
# 5) Compile Third Party software  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Compiling ThirdParty..."  
  setup_environment
  cd $WM_THIRD_PARTY_DIR  
  ./Allwmake -j -l  
  echo "...done."  
  echo ========================================  
fi  
  
# 6) Compile OpenFOAM  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Compiling OpenFOAM..."  
  setup_environment  
  cd $WM_PROJECT_DIR  
  ./Allwmake -j -l  
  echo "...done."  
  echo ========================================  
fi  
  
# 7) Compile FireFOAM  
STEP="$((STEP+1))"
if [ $START_STEP -le $STEP ]; then  
  echo ============ STEP $STEP ================  
  echo "Compiling FireFOAM..."  
  setup_environment  
  cd $INSTALLATION_DIR/$MYFIREFOAM_REPO  
  mkdir build  
  cd build  
  cmake ..  
  make -j"$(nproc)"  
  make install  
  echo "...done."  
  echo ========================================  
fi  
  
echo "Build completed on $(date)"  

