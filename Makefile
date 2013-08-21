# FIXME Parallel processing is not supported in current version
#PARALLEL_FLAGS = -DPARALLEL_SUPPORT -fopenmp
#PARALLEL_LIB = -lgomp
CC = g++
LIBRLCSAPATH = incbwt/
LIBCDSPATH = libcds/
# FIXME -fpermissive is needed for <bcr-demo.o>
CPPFLAGS = -Wall -I$(LIBCDSPATH)includes/ -I$(LIBRLCSAPATH) -g -DMASSIVE_DATA_RLCSA $(PARALLEL_FLAGS) -std=c++0x -fpermissive -O3 -DNDEBUG
LIBCDS = $(LIBCDSPATH)lib/libcds.a
LIBRLCSA = $(LIBRLCSAPATH)/rlcsa.a

INDEXOBJS = CGkArray.o Tools.o HuffWT.o BitRank.o 

all: cgkquery builder

cgkquery: $(LIBCDS) $(LIBRLCSA) $(INDEXOBJS) cgkquery.o
	$(CC) $(CPPFLAGS) -o cgkquery cgkquery.o $(INDEXOBJS) $(LIBCDS) $(LIBRLCSA) $(PARALLEL_LIB)

builder: $(LIBCDS) $(LIBRLCSA) $(INDEXOBJS) builder.o bcr-demo.o
	$(CC) $(CPPFLAGS) -o builder builder.o $(INDEXOBJS) $(LIBCDS)  $(LIBRLCSA) bcr-demo.o

$(LIBCDS):
	@make -C $(LIBCDSPATH)

depend:
	g++ -I$(LIBRLCSAPATH) -I$(LIBCDSPATH)includes/ -std=c++0x -MM *.cpp > dependencies.mk

$(LIBRLCSA):
	@make -C $(LIBRLCSAPATH) library

clean:
	rm -f core *.o *~ builder cgkquery
	@make -C $(LIBCDSPATH) clean
	@make -C $(LIBRLCSAPATH) clean

shallow_clean:
	rm -f core *.o *~ builder cgkquery

include dependencies.mk
