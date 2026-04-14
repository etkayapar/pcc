#! /usr/bin/env bash

FAKECONDAPATH="$HOME/.local/bin/conda"

usage () {
	echo ""
	echo "Usage: create_dummy_conda [auto|</path/to/container/image.simg>]"
	echo "-----"
	echo "If you specify 'auto' as the first argument, then this script will try to find a singularity image"
	echo " inside the directory you run it (should be the pcc directory). If it finds more than one container"
	echo " image then it will exit with an error and will ask you to specify one of the container's path as"
	echo " the first argument to select that specific container." 
	exit 1
}


[ $# -eq 0 ] && usage

if [ $1 == "auto" ]
then
	IFS=$'\n' containers=( $(find . -type f -name '*.simg') )
	n_containers=${#containers[@]}
	if [ $n_containers -gt 1 ]
	then
		echo "ERROR: More than one image found!!!"
		echo "Please identify the correct container image from the found ones below:"
		IFS=$'\n'; printf '%s\n' "${containers[*]}"
		exit 1
	fi

	containerfile=${containers[0]}
else
	containerfile=$1
fi

printf "Selected container path is:  %s\n\n" $containerfile

if ! [ -f $containerfile ]
then
	echo "ERROR: the file $containerfile does not exist!!"
	exit 1
fi

echo "Running 'conda info --json' inside the container to save the JSON string"

apptainer exec --bind="/users,/projappl,/scratch,$TMPDIR,$LOCAL_SCRATCH" ${containerfile} bash -c "conda info --json" | cat <(printf '#! /usr/bin/env bash\n\nJSON=$(cat - <<EOJSON\n') - <(printf '\nEOJSON\n)\n\necho "$JSON"\n') > $FAKECONDAPATH && chmod +x $FAKECONDAPATH

printf "\nFake conda is written to %s .\nMake sure the parent directory (%s) is in your PATH variable. \n\n" $FAKECONDAPATH $(dirname $FAKECONDAPATH) 

echo "Directories that currently are in your path:"
echo "==========================================="
echo $PATH | awk 'BEGIN{RS=":"}{print}'
echo "==========================================="






