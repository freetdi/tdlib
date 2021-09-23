#!/bin/sh

set -e

out=$( mktemp )

./tdecomp --bmd < ${srcdir}/gr/h.gr > $out
./td-validate ${srcdir}/gr/h.gr $out

rm $out
