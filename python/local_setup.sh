#!/bin/bash

printf 'Setup Python2.7 for pyroot'
printf '\n'

export LD_LIBRARY_PATH=/usr/local/Cellar/root/6.12.04_1/lib/root/
export PYTHONPATH=/usr/local/Cellar/root/6.12.04_1/lib/root/
export ROOTSYS=/usr/local/Cellar/root/6.12.04_1/lib/root/

alias python="python2.7"
