#!/usr/bin/env bash

if [[ -n ${1} ]]
then

    if [[ -n ${2} ]]
    then
        exec ${@}
    else
        exec ${1} --help
    fi
else

    echo
    echo "Subread package: high-performance read alignment, quantification and mutation discovery"
    echo "Release: ${SUBREAD_VERSION}"
    echo "http://subread.sourceforge.net/"
    echo
    echo "Usage: singularity run --app subread ${SINGULARITY_CONTAINER} <command>"
    echo
    echo "OR,   provided this alias is defined:"
    echo "      alias subread='singularity run --app subread ${SINGULARITY_CONTAINER}'"
    echo
    echo "Usage: subread <command>"
    echo
    echo "Commands:"
    echo "    Subread"
    echo "    Subjunc"
    echo "    featureCounts"
    echo "    Sublong"
    echo "    exactSNP"
    echo

fi
