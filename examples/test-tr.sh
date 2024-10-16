#!/bin/sh

set -e -x

out=$( mktemp )
alg=--tr

test() {
	echo $alg ${srcdir}/gr/$1.gr
	./tdecomp $alg < ${srcdir}/gr/$1.gr | grep -v incomplete > $out
	./td-validate ${srcdir}/gr/$1.gr $out
}

test h
test path
test path10
test clique3
test clique4

rm $out
