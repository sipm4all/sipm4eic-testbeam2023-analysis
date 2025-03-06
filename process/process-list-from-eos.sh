#! /usr/bin/env bash

MAXDOWNLOAD=1
MAXLIGHTWRITER=16

if [ $# -ne 3 ]; then
    echo " usage: $0 [runlist] [selection] [version] "
    exit 1
fi
RUNLIST=$1
SELECTION=$2
PROCESSVER=$3

PROCESSTAG="process-data-${PROCESSVER}"

mkdir -p log
while read -r run nspills _; do
    [[ ${run} =~ ^#.* ]] && continue
#    if [ -f process/${run}/${PROCESSTAG}/lightdata.root ]; then echo " run ${run} already processed "; continue; fi
    while [ $(ps -ef | grep curl | grep cernbox | wc -l ) -ge ${MAXDOWNLOAD} ]; do sleep 1; done
    while [ $(ps -ef | grep root.exe | grep lightwriter | wc -l ) -ge ${MAXLIGHTWRITER} ]; do sleep 1; done
    echo " --- process-run-from-eos.sh ${run} ${nspills} "
    ./process-run-from-eos.sh ${run} ${nspills} ${SELECTION} ${PROCESSVER} &> log/${run}.process-run-from-eos.log &
    sleep 1
done < ${RUNLIST}

wait
