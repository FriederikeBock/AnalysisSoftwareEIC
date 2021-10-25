#!/usr/bin/env bash

# We need the directory where this is stored so we can determine the other paths,
# regardless of where we call the script.
# NOTE: Using the realpath formulation doesn't seem to work with parsl for whatever
#       reason. I'm sure it could be resolved, but it's not worth messing around with
#       since this alternative seems to work
#currentDir=$(realpath $(dirname "$0"))
# From: https://stackoverflow.com/a/9107028/12907985
if [[ -z "${BASH_SOURCE}" ]]; then
    # zsh
    currentDir="$(dirname $(realpath "$(readlink -f "$0")"))"
else
    # bash
    currentDir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
echo "Loading cpp externals from: ${currentDir}"

# LHAPDF is directly installed in ${currentDir}/install
# Keep the fastjet install separate so we can switch more easily between alibuild fastjet and standalone fastjet
export FASTJET_ROOT="$currentDir/install/fastjet"

export PATH="${currentDir}/install/bin:${FASTJET_ROOT}/bin:${PATH}"
export LD_LIBRARY_PATH="${currentDir}/install/lib:${FASTJET_ROOT}/lib:${LD_LIBRARY_PATH}"
export ROOT_INCLUDE_PATH="${currentDir}/install/include:${FASTJET_ROOT}/include:${ROOT_INCLUDE_PATH}"
export PYTHONPATH="${PYTHONPATH}:${currentDir}/install/lib/python3.9/site-packages"
