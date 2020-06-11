@echo off

echo "Clean windows_visual_studio"

del /Q *.lib

del /Q *.sdf
del /Q /A *.sdf
del /Q /AH *.sdf

del /Q *.opensdf
del /Q /A *.opensdf
del /Q /AH *.opensdf

rmdir /Q /S ".vs"
rmdir /Q /S "Win32"
rmdir /Q /S "x64"

for /D %%d in (*) DO for /D %%e in (%%d\Win32* %%d\x64*) DO rmdir /Q /S "%%e"
for /D %%d in (*) DO for %%e in (%%d\*.exe %%d\*.pdb %%d\*.user) DO del /Q "%%e"
