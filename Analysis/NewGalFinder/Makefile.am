bin_PROGRAMS = gfind

gfind_SOURCES = gfind.c Memory2.c force_spline.mod3.c ost.c \
                subhaloden.mod6.c nrutil.c b2l.c spline.mod2.c \
                nnost.c utils.c mkRtidal.c tsc_omp2.c gsmooth.c

AM_CPPFLAGS = -I$(FFTW)/include
AM_CFLAGS = -g -qopenmp -DINDEX -DVarPM -DXYZDBL -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4 -DDEBUG -DADV -DNMEG=400000L
AM_FCFLAGS = -g -qopenmp -DINDEX -DVarPM -DXYZDBL -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4 -DDEBUG -DADV -mcmodel=medium -DNMEG=400000L

gfind_LDADD = -L$(FFTW)/lib -lsrfftw -lsfftw -lm

CC = mpicc
FC = mpifort
