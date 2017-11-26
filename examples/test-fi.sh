#!/bin/bash

set -e

out=$( mktemp )

./tdecomp --fi < ${srcdir}/gr/h.gr > $out
./td-validate ${srcdir}/gr/h.gr $out

rm $out
