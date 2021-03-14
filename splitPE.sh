#!/bin/bash

SDIR="$( cd "$( dirname "$0" )" && pwd )"

if [ "$VIRTUAL_ENV" != "" ]; then
    VENV_NAME=$(basename $VIRTUAL_ENV)
    if [ "$VENV_NAME" != "scCNATools" ]; then
        . $SDIR/venv/scCNATools/bin/activate
    fi
else
    . $SDIR/venv/scCNATools/bin/activate
fi

$SDIR/SplitPE/splitPE.py $@
