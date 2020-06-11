#!/bin/sh

rm -f $1.log
rm -f $1.err
rm -f $1.post.msh
rm -f $1.???.post.msh
rm -f $1.geometry.dat

# OutputFile: $1.log
# ErrorFile: $1.err

mv $1.dat $1.geometry.dat

$3/Coloring $1 2 $1.log 2> $1.err
