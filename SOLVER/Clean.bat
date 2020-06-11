@echo off

pushd build\GCC
call Clean.bat
popd

pushd build\ICC
call Clean.bat
popd

pushd "build\Visual Studio"
call Clean.bat
popd

pushd "build\Visual Studio 2015"
call Clean.bat
popd

pushd gid\examples
for /D %%d in (*.gid) DO for %%f in (%%d\*.dat %%d\*.res %%d\*.log %%d\*.err %%d\*.msh %%d\*.vv %%d\*.png %%d\*.rdr %%d\*.bin %%d\*post.mat) DO del /Q %%f
popd

pushd gid\problemtypes
for /D %%d in (*.gid) DO for %%f in (%%d\*.exe %%d\*.pdb) DO del /Q %%f
popd

pushd tools
del /Q *.exe
del /Q *.pdb
del /Q *.mat
popd
