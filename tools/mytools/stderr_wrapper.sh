TMPFILE=$(mktemp) || exit 1

"$@" 2> $TMPFILE

EXITCODE=$?
# Exitcode != 0 ?
if [ "$EXITCODE" -ne "0" ]; then
       echo "Error code $EXITCODE" >&2
       cat -n -v $TMPFILE >&2
fi
rm $TMPFILE

exit $EXITCODE
