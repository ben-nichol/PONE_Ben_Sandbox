_PONE_SRC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
_PYTHONPATH=$_PONE_SRC:$PYTHONPATH

export PYTHONPATH=$_PYTHONPATH 
export PONESRCDIR=$_PONE_SRC 

printf "PONE environment has:\n"
printf "   PONE_SRC     = %s\n" $PONESRCDIR

