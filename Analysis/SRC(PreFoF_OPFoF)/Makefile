#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  opfof.exe 

#--- Objects to build

OBJS		= opfof.o Time.o Treewalk.fof.ordered.o  read.o  

			# fda4_new.o fda4_split.o tsc_c.o

#			fmod.o fda4_new.o fda4_split.o tsc_c.o \
#indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#hypersort.o\
#slicelightcone.final.o
#			kjhtpm.final.o indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#			findindex.o maxextent.o glacial.o survey.o surveywrite2.o
#			memory.o	tpm.o domdecomp.o  kjhgetvisitors.o \
# 			tree.o  quickbin.o group.o groupgrav.o  \
#			kjhtpm.o Treewalk.mod6.o force_spline.mod2.o
#	hypersort.o \

OS = LINUX_GCC



#--- Targets

CAMBDIR = 
#CAMBDIR = CAMB
SUBDIRS = $(CAMBDIR)

all: 
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) all);\
	done 
	$(MAKE) this

this: $(OBJS)
	rm -f $(PROGS) 
	$(MAKE) $(PROGS)
new: 
	$(MAKE) clean
	$(MAKE) all

clean: 
	rm -f $(OBJS)
	rm -f $(OBJPOW)
	rm -f $(OBJVIEW)
	rm -f $(PROGS)
	rm -f $(MKPOW)
	rm -f $(VIEWPOW)
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) clean);\
	done 


#mkpower:
#	rm -f $(MKPOW) 
#	rm -f $(OBJPOW) 
#	$(MAKE) $(MKPOW)

default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
