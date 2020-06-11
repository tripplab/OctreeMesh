#!/bin/sh

rm -f $1.log
rm -f $1.err
rm -f $1.wrn
rm -f $1.post.res
rm -f $1.post.msh
rm -f $1.???.post.res
rm -f $1.???.post.msh
rm -f $1.geometry.dat
rm -f $1.problem.dat
rm -f $1.solver.dat

# OutputFile: $1.log
# ErrorFile: $1.err
# WarningFile: $1.wrn

mv $1.dat $1.geometry.dat
mv $1-1.dat $1.problem.dat
mv $1-2.dat $1.solver.dat

$3/Solid $1 2 $1.log 2> $1.err
