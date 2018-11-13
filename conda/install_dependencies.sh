#!/bin/bash
# Stop on error
set -e

CONDA_ENV=circAlign

SH_SCRIPT_DIR=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

REQ_TXT=${SH_SCRIPT_DIR}/requirements.txt

if which conda; then
  echo "=== Found Conda ($(conda --version))."
else
  echo "=== Conda does not exist on your system. Please install Conda first."
  echo "https://conda.io/docs/user-guide/install/index.html#regular-installation"
  exit 1
fi

if conda env list | grep -wq ${CONDA_ENV}; then
  echo "=== Pipeline's circAlign Conda env (${CONDA_ENV}) already exists."
  echo "=== Please remove it first (uninstall_dependencies.sh)"
  exit 2
fi

echo "=== Installing packages for CircAlign env..."
conda create -n ${CONDA_ENV} --file ${REQ_TXT} -y -c bioconda -c conda-forge -c defaults -c r

if [[ -f ~/.bashrc ]]; then
    echo "export PYTHONPATH=${SH_SCRIPT_DIR}/.." >> "~/.bashrc"
elif [[ -f ~/.bash_profile ]]; then
    echo "export PYTHONPATH=${SH_SCRIPT_DIR}/.." >> "~/.bash_profile"
else
    echo "Could not add circAlign to bash profile, do so manually to avoid issues"
    export PYTHONPATH="${SH_SCRIPT_DIR}/.."
fi

echo "=== All done."

