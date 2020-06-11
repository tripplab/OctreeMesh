.SILENT:
.PHONY: all clean

SRCDIR:=../../../source

CPPFLAGS:=-I$(SRCDIR)/METIS
ARFLAGS:=rv

#CFLAGS:=-g
CFLAGS:=-O2 -DNDEBUG
CC:=gcc

all: ../libMETIS.a

clean:
	echo CLEAN
	rm -f *.o ../libMETIS.a

######################################################################
# METIS ##############################################################
######################################################################
../libMETIS.a: balance.o bucketsort.o ccgraph.o coarsen.o compress.o debug.o estmem.o fm.o fortran.o frename.o graph.o initpart.o kmetis.o kvmetis.o kwayfm.o kwayrefine.o kwayvolfm.o kwayvolrefine.o match.o mbalance.o mbalance2.o mcoarsen.o memory.o mesh.o meshpart.o mfm.o mfm2.o mincover.o minitpart.o minitpart2.o mkmetis.o mkwayfmh.o mkwayrefine.o mmatch.o mmd.o mpmetis.o mrefine.o mrefine2.o mutil.o myqsort.o ometis.o parmetis.o pmetis.o pqueue.o refine.o separator.o sfm.o srefine.o stat.o subdomains.o timing.o util.o
	echo AR $@
	$(AR) rcs $@ $+

balance.o: $(SRCDIR)/METIS/balance.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

bucketsort.o: $(SRCDIR)/METIS/bucketsort.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

ccgraph.o: $(SRCDIR)/METIS/ccgraph.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

coarsen.o: $(SRCDIR)/METIS/coarsen.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

compress.o: $(SRCDIR)/METIS/compress.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

debug.o: $(SRCDIR)/METIS/debug.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

estmem.o: $(SRCDIR)/METIS/estmem.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

fm.o: $(SRCDIR)/METIS/fm.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

fortran.o: $(SRCDIR)/METIS/fortran.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

frename.o: $(SRCDIR)/METIS/frename.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

graph.o: $(SRCDIR)/METIS/graph.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

initpart.o: $(SRCDIR)/METIS/initpart.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kmetis.o: $(SRCDIR)/METIS/kmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kvmetis.o: $(SRCDIR)/METIS/kvmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kwayfm.o: $(SRCDIR)/METIS/kwayfm.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kwayrefine.o: $(SRCDIR)/METIS/kwayrefine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kwayvolfm.o: $(SRCDIR)/METIS/kwayvolfm.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

kwayvolrefine.o: $(SRCDIR)/METIS/kwayvolrefine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

match.o: $(SRCDIR)/METIS/match.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mbalance.o: $(SRCDIR)/METIS/mbalance.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mbalance2.o: $(SRCDIR)/METIS/mbalance2.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mcoarsen.o: $(SRCDIR)/METIS/mcoarsen.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

memory.o: $(SRCDIR)/METIS/memory.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mesh.o: $(SRCDIR)/METIS/mesh.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

meshpart.o: $(SRCDIR)/METIS/meshpart.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mfm.o: $(SRCDIR)/METIS/mfm.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mfm2.o: $(SRCDIR)/METIS/mfm2.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mincover.o: $(SRCDIR)/METIS/mincover.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

minitpart.o: $(SRCDIR)/METIS/minitpart.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

minitpart2.o: $(SRCDIR)/METIS/minitpart2.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mkmetis.o: $(SRCDIR)/METIS/mkmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mkwayfmh.o: $(SRCDIR)/METIS/mkwayfmh.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mkwayrefine.o: $(SRCDIR)/METIS/mkwayrefine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mmatch.o: $(SRCDIR)/METIS/mmatch.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mmd.o: $(SRCDIR)/METIS/mmd.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mpmetis.o: $(SRCDIR)/METIS/mpmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mrefine.o: $(SRCDIR)/METIS/mrefine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mrefine2.o: $(SRCDIR)/METIS/mrefine2.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

mutil.o: $(SRCDIR)/METIS/mutil.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

myqsort.o: $(SRCDIR)/METIS/myqsort.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

ometis.o: $(SRCDIR)/METIS/ometis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

parmetis.o: $(SRCDIR)/METIS/parmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

pmetis.o: $(SRCDIR)/METIS/pmetis.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

pqueue.o: $(SRCDIR)/METIS/pqueue.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

refine.o: $(SRCDIR)/METIS/refine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

separator.o: $(SRCDIR)/METIS/separator.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

sfm.o: $(SRCDIR)/METIS/sfm.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

srefine.o: $(SRCDIR)/METIS/srefine.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

stat.o: $(SRCDIR)/METIS/stat.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

subdomains.o: $(SRCDIR)/METIS/subdomains.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

timing.o: $(SRCDIR)/METIS/timing.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

util.o: $(SRCDIR)/METIS/util.c
	echo CC $<
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

