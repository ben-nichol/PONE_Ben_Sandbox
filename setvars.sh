#!/bin/sh

_PONE_SRC=$(cd -- "$(dirname -- $0)" && pwd)
_PYTHONPATH="${PYTHONPATH}:${_PONE_SRC}"
_LD_LIBRARY_PATH=/usr/local/OpenBLAS:$LD_LIBRARY_PATH

export PONESRCDIR="$_PONE_SRC"
export PYTHONPATH="$_PYTHONPATH"

#Check if pone_offline version matches I3 version
I3VERSION=""
RED='\033[1;31m'
GREEN='\033[1;32m'
NC='\033[0m'
I3_CMAKECACHE="$I3_SRC/icetray/CMakeCache.txt"
if [ -f "$I3_CMAKECACHE" ]; then
	I3VERSION=$(grep '^CMAKE_PROJECT_VERSION:STATIC=' "$I3_CMAKECACHE" | cut -d'=' -f2)
else
	echo "${RED}WARNING: CMakeCache.txt not found at $I3_CMAKECACHE${NC}" 1>&2
fi
PTAG=""
if [ -f "${_PONE_SRC}/.versiontag" ]; then
	PTAG=$(cat "${_PONE_SRC}/.versiontag")
fi
if [ -n "$I3VERSION" ] && [ -n "$PTAG" ] && [ "$I3VERSION" != "$PTAG" ]; then
	echo "${RED}WARNING: PONE_OFFLINE EXPECTS ICETRAY VERSION $PTAG, but found $I3VERSION${NC}" 1>&2
elif [ -z "$I3VERSION" ] || [ -z "$PTAG" ]; then
	echo "${RED} WARNING: MISSING VERSION"
	echo "${RED} PONE_OFFLINE: ${PTAG}"
	echo "${RED} ICETRAY: ${I3VERSION}${NC}"
else
	echo "${GREEN}PONE_OFFLINE FOUND EXPECTED ICETRAY VERSION${NC}" 1>&2
fi


if [ "$#" -eq 1 ]; then
	echo "P-ONE Offline environment has:"
	echo "   I3_SRC     = $I3_SRC"
	echo "   I3_BUILD   = $I3_BUILD"
	echo "   PONESRCDIR = $PONESRCDIR"
fi

exec "$@"
