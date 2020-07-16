#!/bin/bash
# if errors emerge (usually relating to allensdk) execute the following:
# xcode-select --install # MacOS Xcode developer tools
# or installing homebrew should automatically install xcode conmmand line dev tools
conda env create -f environment.yml # or change to environment-big.yml
conda activate msbrainpy
pip install allensdk
pip install requests
python setup.py install
