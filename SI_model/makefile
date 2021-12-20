#
# This is the makefile for libcash2-s.so
#

###############################################
# Specify where you want to install the library
INSTALL_DIRECTORY = $(HOME)/local

##############################
# Specify which program to use
CC = gcc
SED = sed
INSTALL = install
RM = /bin/rm -f
RMDIR = rmdir
#########################################################
# Do not modify below unless you know what you want to do
OBJCASH = arithmetic.o basic.o color.o filter.o io.o logical.o \
margolus.o movie.o neighbors.o noise.o png.o ps.o random.o shift.o\
x11.o
OBJCASH2 = cash2.o mersenne.o 
OBJMAIN = cash2-s.o

CCOPT = -fPIC -c -O3 -Wall
#CCOPT = -fPIC -c -ggdb -Wall

LDIR = -L/usr/X11R6/lib
LIBS = -lpng -lz -lX11 -lm -lgrace_np

all: libcash2-s.a make.makefile

libcash2-s.so: $(OBJCASH) $(OBJCASH2) $(OBJMAIN)
	$(CC) -shared -Wl,-soname,libcash2-s.so -o libcash2-s.so $(OBJCASH) $(OBJCASH2) $(OBJMAIN) -lc -lpng -lz -lX11 -lm -lgrace_np -L/usr/X11R6/lib

libcash2-s.a: $(OBJCASH) $(OBJCASH2) $(OBJMAIN)
	rm -f libcash2-s.a
	ar svr libcash2-s.a $(OBJCASH) $(OBJCASH2) $(OBJMAIN)

make.makefile:
	$(SED) -e 's|$$(HOME)|$(INSTALL_DIRECTORY)|g' makefile.project > Model/makefile

%.o: %.c
	$(CC) $(CCOPT) $<

$(OBJCASH): cash2003.h  makefile
$(OBJCASH2): cash2.h mersenne.h  makefile
$(OBJMAIN): cash2003.h cash2.h mersenne.h cash2-s.h makefile

install: all
	$(INSTALL) -d $(INSTALL_DIRECTORY)
	$(INSTALL) -d $(INSTALL_DIRECTORY)/include $(INSTALL_DIRECTORY)/lib
	$(INSTALL) -m 644 cash2003.h cash2.h cash2-s.h mersenne.h $(INSTALL_DIRECTORY)/include
	$(INSTALL) libcash2-s.a $(INSTALL_DIRECTORY)/lib

uninstall:
	$(RM) $(INSTALL_DIRECTORY)/include/cash2003.h
	$(RM) $(INSTALL_DIRECTORY)/include/cash2.h
	$(RM) $(INSTALL_DIRECTORY)/include/cash2-s.h
	$(RM) $(INSTALL_DIRECTORY)/include/mersenne.h
	$(RMDIR) --ignore-fail-on-non-empty $(INSTALL_DIRECTORY)/include
	$(RM) $(INSTALL_DIRECTORY)/lib/libcash2-s.so
	$(RM) $(INSTALL_DIRECTORY)/lib/libcash2-s.a
	$(RMDIR) --ignore-fail-on-non-empty $(INSTALL_DIRECTORY)/lib
	$(RMDIR) --ignore-fail-on-non-empty $(INSTALL_DIRECTORY)

clean:
	$(RM) *.o libcash2-s.so libcash2-s.a
