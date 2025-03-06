#! /usr/bin/env bash

### standard for version / selection (to be defined)
### v0.0 - no selection
### v0.1 - trigger 0 selection
### v0.2 - timing selection

CERNBOXDIR="https://cernbox.cern.ch/remote.php/dav/public-files/jkHfV7ISSRVUbip/physics/all"
LOCALDIR="/home/eic/DATA/2024-testbeam/process"

if [ $# -ne 4 ]; then
    echo " usage: $0 [runname] [nspills] [selection] [version] "
    exit 1
fi
RUNNAME=$1
NSPILLS=$2
SELECTION=$3
PROCESSVER=$4

main() {
    RUNNAME=$1
    NSPILLS=$2
    mkdir -p ${LOCALDIR}

    echo " --- moving to ${LOCALDIR} "
    cd ${LOCALDIR}
    download_data ${RUNNAME}
    process_data ${RUNNAME} ${NSPILLS}
    cd - &> /dev/null
}

download_data() {
    RUNNAME=$1

    echo " --- downloading ${RUNNAME} data from cernbox "
    curl ${CERNBOXDIR}/${RUNNAME}.tgz --output ${RUNNAME}.tgz || exit 1
    echo " --- uncompress ${RUNNAME} data "
    tar zxvf ${RUNNAME}.tgz &> /dev/null || exit 1
    echo " --- remove ${RUNNAME} tarball "
    rm -f ${RUNNAME}.tgz || exit 1
    
    cd ${RUNNAME}
    for KC in kc705*; do
	echo " --- uncompress ${RUNNAME} ${KC} decoded data "
	cd ${KC}
#	tar zxvf raw.tgz &> /dev/null
	tar zxvf decoded.tgz &> /dev/null
	rm -f *.tgz
	cd - &> /dev/null
    done
    cd - &> /dev/null
}

process_data() {
    RUNNAME=$1
    NSPILLS=$2
    
    RUNDIR=${LOCALDIR}/${RUNNAME}
    SOFTDIR="/home/eic/alcor/alcor-utils/pdu/measure/sipm4eic-testbeam2023-analysis"
    PROCESSTAG="process-data-${PROCESSVER}"
    PROCESSDIR="${RUNDIR}/${PROCESSTAG}"
    mkdir -p ${PROCESSDIR}

    ### create lightdata
    
    LIGHTWRITER="${SOFTDIR}/macros/lightwriter.C"
    LIGHTDATA="${PROCESSDIR}/lightdata.root"
    FINEDATA="${PROCESSDIR}/finedata.root"

    echo " --- create ${RUNNAME} lightdata "
    echo " --- process lightwriter "
    root -b -l <<-EOF || exit 1
        .L ${LIGHTWRITER}
	apply_selection["${SELECTION}"] = true;    
    	lightwriter("${RUNDIR}", "${LIGHTDATA}", "${FINEDATA}", ${NSPILLS})    
EOF
    
    ### create rolloverdata
    
    ROLLOVERWRITER="${SOFTDIR}/macros/rolloverwriter.C"
    ROLLOVERDATA="${PROCESSDIR}/rolloverdata.root"

    echo " --- create ${RUNNAME} rolloverdata "
    root -b -q -l "${ROLLOVERWRITER}(\"${RUNDIR}\", \"${ROLLOVERDATA}\")" || exit 1

    ### create recodata
    
    RECOWRITER="${SOFTDIR}/macros/recowriter.C"
    RECODATA="${PROCESSDIR}/recodata.root"

    echo " --- create ${RUNNAME} recodata "
    root -b -q -l "${RECOWRITER}(\"${LIGHTDATA}\", \"${RECODATA}\")" || exit 1

    ### delete everything but process dir

    echo " --- delete ${RUNNAME} decoded data "
    cd ${RUNDIR}
    for DIR in *; do
	if [[ ${DIR} == process-data* ]]; then continue; fi
	rm -rf ${DIR}
    done

    cd - &> /dev/null
}

main ${RUNNAME} ${NSPILLS}
