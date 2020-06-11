@echo off

del /Q %1.log
del /Q %1.err
del /Q %1.wrn
del /Q %1.post.res
del /Q %1.post.msh
del /Q %1.???.post.res
del /Q %1.???.post.msh
del /Q %1.geometry.dat
del /Q %1.problem.dat
del /Q %1.solver.dat

rem OutputFile: %1.log
rem ErrorFile: %1.err
rem WarningFile: %1.wrn

ren %1.dat %1.geometry.dat
ren %1-1.dat %1.problem.dat
ren %1-2.dat %1.solver.dat

%3\Solid %1 2 %1.log 2> %1.err
