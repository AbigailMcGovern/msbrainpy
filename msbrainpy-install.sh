#!/bin/bash
# if errors emerge (usually relating to allensdk) execute the following:
# xcode-select --install # MacOS Xcode developer tools
conda env create -f environment.yml # or change to environment-big.yml
conda activate msbrainpy
pip install allensdk
pip install requests
sudo python setup.py install
