@echo off

del /Q %1.log
del /Q %1.err
del /Q %1.post.res
del /Q %1.post.msh
del /Q %1.???.post.res
del /Q %1.???.post.msh
del /Q %1.geometry.dat
del /Q %1.problem.dat

rem OutputFile: %1.log
rem ErrorFile: %1.err

ren %1.dat %1.geometry.dat
ren %1-1.dat %1.problem.dat

%3\Partitioning %1 2 %1.log 2> %1.err
