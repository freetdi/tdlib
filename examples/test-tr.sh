#!/bin/sh

set -e

out=$( mktemp )

alg=--tr

./tdecomp $alg < ${srcdir}/gr/h.gr > $out
./td-validate ${srcdir}/gr/h.gr $out
./tdecomp $alg < ${srcdir}/gr/path.gr | grep -v incomplete > $out
./td-validate ${srcdir}/gr/path.gr $out
./tdecomp $alg < ${srcdir}/gr/path10.gr | grep -v incomplete > $out
./td-validate ${srcdir}/gr/path10.gr $out
./tdecomp $alg < ${srcdir}/gr/clique3.gr | grep -v incomplete > $out
./td-validate ${srcdir}/gr/clique3.gr $out
./tdecomp $alg < ${srcdir}/gr/clique4.gr | grep -v incomplete > $out
./td-validate ${srcdir}/gr/clique4.gr $out

rm $out
