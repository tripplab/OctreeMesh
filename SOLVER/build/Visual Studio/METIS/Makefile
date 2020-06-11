!INCLUDE ..\BuildFlags.mk

SRCDIR=..\..\..\source

!IF $(USE_DEBUG)==1
CPPFLAGS=/I$(SRCDIR)\METIS
CFLAGS=/nologo /Zi
!ELSE
CPPFLAGS=/I$(SRCDIR)\METIS /DNDEBUG
CFLAGS=/nologo /O2
!ENDIF

all: ..\METIS.lib

clean:
	del /Q *.obj *.pdb ..\METIS.lib

######################################################################
# METIS ##############################################################
######################################################################
..\METIS.lib: balance.obj bucketsort.obj ccgraph.obj coarsen.obj compress.obj debug.obj estmem.obj fm.obj fortran.obj frename.obj graph.obj initpart.obj kmetis.obj kvmetis.obj kwayfm.obj kwayrefine.obj kwayvolfm.obj kwayvolrefine.obj match.obj mbalance.obj mbalance2.obj mcoarsen.obj memory.obj mesh.obj meshpart.obj mfm.obj mfm2.obj mincover.obj minitpart.obj minitpart2.obj mkmetis.obj mkwayfmh.obj mkwayrefine.obj mmatch.obj mmd.obj mpmetis.obj mrefine.obj mrefine2.obj mutil.obj myqsort.obj ometis.obj parmetis.obj pmetis.obj pqueue.obj refine.obj separator.obj sfm.obj srefine.obj stat.obj subdomains.obj timing.obj util.obj
	$(AR) /NOLOGO /OUT:$@ $**

balance.obj: $(SRCDIR)\METIS\balance.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

bucketsort.obj: $(SRCDIR)\METIS\bucketsort.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

ccgraph.obj: $(SRCDIR)\METIS\ccgraph.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

coarsen.obj: $(SRCDIR)\METIS\coarsen.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

compress.obj: $(SRCDIR)\METIS\compress.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

debug.obj: $(SRCDIR)\METIS\debug.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

estmem.obj: $(SRCDIR)\METIS\estmem.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

fm.obj: $(SRCDIR)\METIS\fm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

fortran.obj: $(SRCDIR)\METIS\fortran.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

frename.obj: $(SRCDIR)\METIS\frename.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

graph.obj: $(SRCDIR)\METIS\graph.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

initpart.obj: $(SRCDIR)\METIS\initpart.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kmetis.obj: $(SRCDIR)\METIS\kmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kvmetis.obj: $(SRCDIR)\METIS\kvmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kwayfm.obj: $(SRCDIR)\METIS\kwayfm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kwayrefine.obj: $(SRCDIR)\METIS\kwayrefine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kwayvolfm.obj: $(SRCDIR)\METIS\kwayvolfm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

kwayvolrefine.obj: $(SRCDIR)\METIS\kwayvolrefine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

match.obj: $(SRCDIR)\METIS\match.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mbalance.obj: $(SRCDIR)\METIS\mbalance.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mbalance2.obj: $(SRCDIR)\METIS\mbalance2.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mcoarsen.obj: $(SRCDIR)\METIS\mcoarsen.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

memory.obj: $(SRCDIR)\METIS\memory.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mesh.obj: $(SRCDIR)\METIS\mesh.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

meshpart.obj: $(SRCDIR)\METIS\meshpart.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mfm.obj: $(SRCDIR)\METIS\mfm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mfm2.obj: $(SRCDIR)\METIS\mfm2.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mincover.obj: $(SRCDIR)\METIS\mincover.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

minitpart.obj: $(SRCDIR)\METIS\minitpart.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

minitpart2.obj: $(SRCDIR)\METIS\minitpart2.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mkmetis.obj: $(SRCDIR)\METIS\mkmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mkwayfmh.obj: $(SRCDIR)\METIS\mkwayfmh.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mkwayrefine.obj: $(SRCDIR)\METIS\mkwayrefine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mmatch.obj: $(SRCDIR)\METIS\mmatch.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mmd.obj: $(SRCDIR)\METIS\mmd.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mpmetis.obj: $(SRCDIR)\METIS\mpmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mrefine.obj: $(SRCDIR)\METIS\mrefine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mrefine2.obj: $(SRCDIR)\METIS\mrefine2.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

mutil.obj: $(SRCDIR)\METIS\mutil.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

myqsort.obj: $(SRCDIR)\METIS\myqsort.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

ometis.obj: $(SRCDIR)\METIS\ometis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

parmetis.obj: $(SRCDIR)\METIS\parmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

pmetis.obj: $(SRCDIR)\METIS\pmetis.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

pqueue.obj: $(SRCDIR)\METIS\pqueue.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

refine.obj: $(SRCDIR)\METIS\refine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

separator.obj: $(SRCDIR)\METIS\separator.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

sfm.obj: $(SRCDIR)\METIS\sfm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

srefine.obj: $(SRCDIR)\METIS\srefine.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

stat.obj: $(SRCDIR)\METIS\stat.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

subdomains.obj: $(SRCDIR)\METIS\subdomains.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

timing.obj: $(SRCDIR)\METIS\timing.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

util.obj: $(SRCDIR)\METIS\util.c
	$(CC) $(CPPFLAGS) $(CFLAGS) /Fo$@ /c $?

