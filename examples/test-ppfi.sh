#!/bin/sh

set -e

alg=--ppfi

out=$( mktemp )

./tdecomp $alg < ${srcdir}/gr/h.gr | grep -v incomplete > $out
./td-validate ${srcdir}/gr/h.gr $out

rm $out
