#!/bin/sh

# There are two cases we want to handle here: 
# 1. The IceTray environment is not loaded, so we should try to run it, injecting our variables
# 2. We are already inside a loaded IceTray environment, and only need to tweak a couple of variables
#
# In either case, we want to preserve the ability to run a command specified as arguments.

if [ "_$PONESRCDIR" != "_" ]; then
	echo "*************************************************************" 1>&2
	echo "The pone_offline environment appears to already be configured" 1>&2
	echo "*************************************************************" 1>&2
	exit 1
fi

# At the moment, pone_offline only runs 'in place', without any installation step.
# As a result, wherever this script is is its root.

_PONE_SRC=$(cd -- "$(dirname -- $0)" && pwd)

if [ -z "$1" ]; then
	NEW_SHELL=$SHELL
else
	NEW_SHELL=$1
	shift
fi

# The IceTray environment sets I3_SHELL to identify itself
if [ ! "$I3_SHELL" ]; then
	# environment not loaded, so we will attempt to do it
	if [ "_$(command -v "icetray-shell")" = "_" ]; then
		echo "****************************************************************" 1>&2
		echo "icetray-shell not found in \$PATH, unable to start automatically" 1>&2
		echo "****************************************************************" 1>&2
		exit 2
	fi
	exec icetray-shell "${_PONE_SRC}/setvars.sh" "$NEW_SHELL" "$@"
else
	exec "${_PONE_SRC}/setvars.sh" "$NEW_SHELL" "$@"

fi
