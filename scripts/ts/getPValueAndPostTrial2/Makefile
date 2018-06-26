 ######################################################################
 # Project: Examples of how to use TMVA                               #
 ###################################################################### 

prefix=$(HOME)/ADC/standard_analysis/

exec_prefix=$(prefix)
MVAINSTALL=$(prefix)

MAKEFLAGS = --no-print-directory -r -s

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)


INCLUDE = $(ROOTCFLAGS) -I /usr/local_gtb/gmp/include/
LIBS    = $(ROOTLIBS) -L/usr/local_gtb/gmp/lib/ -lgmp

BIN1 =  getPValueAndPostTrial
#BINS2 = BKG_TMVAINPUT

CCFLAGS   = -O0 -m128bit-long-double
AR       = ar cqs
TAR      = tar -cf
GZIP     = gzip -9f
COPY     = cp -f
COPY_FILE= $(COPY) -p
COPY_DIR = $(COPY) -pR
DEL_FILE = rm -f
SYMLINK  = ln -sf
DEL_DIR  = rm -rf
MOVE     = mv -f
CHK_DIR_EXISTS= test -d
MKDIR    = mkdir -p
CVS = cvs 
LINK     = g++
LFLAGS = -shared -Wl,-soname,$(TARGET0) -Wl,-rpath,$(DESTDIR)

all: $(BIN1) $(BINS2)

$(BIN1): % : %.cpp 
	@echo -n "Building $@ ... "
	@echo	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -o $@
	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -o $@
	@echo "Done"

$(BINS2): % : %.cpp 
	@echo -n "Building $@ ... "
	@echo	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -o $@
	$(CXX) $(CCFLAGS) $< $(INCLUDE) $(LIBS) -o $@
	@echo "Done"

clean:
	rm -f $(BIN1)
	rm -f $(BIN2)

install: all
	$(COPY_FILE) ./$(BINS) $(MVAINSTALL)/bin

install_test: all
	$(COPY_FILE) ./$(BINS2) ~/TMVA
	$(COPY_FILE) ./fileinput.conf ~/TMVA
	
