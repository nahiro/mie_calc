#!/bin/bash

if [ "${#}" != "5" ] ; then
  echo "Usage: auto_submit.sh DATDIR CALDIR SMOD IATM IHAZ" >&2
  exit
fi
DATDIR=${1}
CALDIR=${2}
SMOD=${3}
IATM=${4}
IHAZ=${5}
INAM=input.dat
CNAM=conf.dat

cd ${DATDIR}
if [ ! -f ${CNAM} ] ; then
  echo "No such file >>> ${CNAM}" >&2
  exit
fi
if [ ! -f ${INAM} ] ; then
  echo "No such file >>> ${INAM}" >&2
  exit
fi

nohup aeros_mix -ccc -vv -C ${CNAM} -D ${CALDIR} -s ${SMOD} -A ${IATM} -a ${IHAZ} <${INAM} >>nohup.out 2>&1 &
