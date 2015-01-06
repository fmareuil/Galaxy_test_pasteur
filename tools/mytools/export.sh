#!/bin/bash

while getopts "i:f:n:e:o:" optionName; do
case "$optionName" in

i) INPUT="$OPTARG";;
f) INPUT_NAME="$OPTARG";;
n) NAME="$OPTARG";;
e) EXPORT_DIR="$OPTARG";;
o) OUTPUT="$OPTARG";;
esac
done

EXPORT_DIR=$EXPORT_DIR/../projets/export
a=`echo $NAME | sed 's/[A-Za-z0-9_-]//g' | sed 's/\.//'`
if [ -z "${a}" ]; then 
    if [ -e "$EXPORT_DIR/$NAME" ]; then
        echo "$EXPORT_DIR/$NAME already exists: delete file or change the name" >&2
	exit 1
    elif [ `du -c $EXPORT_DIR | grep 'total' | awk '{print $1}'` -gt 104857600 ]; then
        echo "$EXPORT_DIR is too large: move or delete files \nif you do not have permission,please contact an admin" >&2
	exit 1
    else
        echo -e "command : cp $INPUT_NAME $EXPORT_DIR/$NAME\n" > $OUTPUT
        cp "$INPUT" "$EXPORT_DIR/$NAME"
        chmod 660 "$EXPORT_DIR/$NAME"
        echo -e "$INPUT_NAME \nhave been copied in \n$EXPORT_DIR \nand named \n$NAME \nWARNING: this tool duplicates the data, remember to delete redundant data.">> $OUTPUT
    fi
else
    e=0
    while [ $e -lt `echo ${#a}` ]; do 
        if [ ${a:$e:$e+1} == "." ]; then
            echo "too many points ('.' character) in the filename" >&2
            exit 1
        else
            echo "wrong character in the filename" >&2
            exit 1
	fi ;
	e=$(($e+1))
    done
fi
