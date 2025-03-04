#!/bin/sh

set -e -x

out=$( mktemp )
alg=--ta

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
test sage_bug_39404

rm $out
