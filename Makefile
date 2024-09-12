CXX=g++
ROOT_FLAGS=`root-config --cflags --libs `
G4_FLAGS=`geant4-config --libs`

GSL_FLAGS=-lgsl

OXSX_ROOT = $(OXO_DIR)
OXSX_INC=$(OXSX_ROOT)/include
OXSX_LIB_DIR=$(OXSX_ROOT)/build 
OXSX_LIB_NAME=oxsx

H5_LIBS = -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5

SRC_FILES := $(wildcard src/*.cc)
OBJ_FILES := $(addprefix build/, $(notdir $(SRC_FILES:.cc=.o)))

PREFIX ?= /usr/local/bin

INC_DIR=src
LIB_DIR=lib
LIB_NAME=antinullh

LIB=$(LIB_DIR)/lib$(LIB_NAME).a

all: bin/make_trees bin/llh_scan #bin/make_plots bin/make_pdfs bin/fit_dataset bin/build_asimov bin/make_plots bin/llh_scan bin/auto_corrs

bin/fit_dataset: fit_dataset.cc $(LIB)
	mkdir -p bin
	$(CXX)  fit_dataset.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) ${GSL_FLAGS} -larmadillo -lMinuit2 -o $@

bin/make_plots: make_plots.cc $(LIB)
	mkdir -p bin
	$(CXX)  make_plots.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

bin/up_count_lim: up_count_lim.cc $(LIB)
	mkdir -p bin
	$(CXX)  up_count_lim.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) $(GSL_FLAGS) -larmadillo -o $@

bin/make_pdfs: make_pdfs.cc $(LIB)
	mkdir -p bin
	$(CXX)  make_pdfs.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

bin/make_trees: make_trees.cc $(LIB)
	mkdir -p bin
	$(CXX)  make_trees.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) -larmadillo -o $@

bin/build_asimov: build_asimov.cc $(LIB)
	mkdir -p bin
	$(CXX)  build_asimov.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) ${GSL_FLAGS} -larmadillo -o $@

bin/auto_corrs: auto_corrs.cc $(LIB)
	mkdir -p bin
	$(CXX)  auto_corrs.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) ${GSL_FLAGS} -larmadillo -o $@

bin/llh_scan: llh_scan.cc $(LIB)
	mkdir -p bin
	$(CXX)  llh_scan.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) ${GSL_FLAGS} -larmadillo -o $@

$(LIB) : $(OBJ_FILES)
	mkdir -p $(LIB_DIR)
	ar rcs  $@ $^

build/%.o : src/%.cc
	mkdir -p build
	$(CXX) -c -w $< -I$(OXSX_INC) -Isrc/ -w $(ROOT_FLAGS) $(G4_FLAGS) -o $@

install:
	ln -sf `readlink -f bin/make_pdfs` $(PREFIX)
	ln -sf `readlink -f bin/make_trees` $(PREFIX)
	ln -sf `readlink -f bin/make_plots` $(PREFIX)
	ln -sf `readlink -f bin/fit_dataset` $(PREFIX)
	ln -sf `readlink -f bin/fit_dataset_batch` $(PREFIX)
	ln -sf `readlink -f bin/build_asimov` $(PREFIX)
	ln -sf `readlink -f bin/auto_corrs` $(PREFIX)
	ln -sf `readlink -f bin/llh_scan` $(PREFIX)
	chmod +x bin/make_pdfs
	chmod +x bin/make_plots
	chmod +x bin/make_trees
	chmod +x bin/fit_dataset
	chmod +x bin/fit_dataset_batch
	chmod +x bin/build_asimov
	chmod +x bin/auto_corrs
	chmod +x bin/llh_scann

clean:
	rm -f bin/make_pdfs
	rm -f bin/make_plots
	rm -f bin/make_trees
	rm -f bin/fit_dataset
	rm -f bin/build_asimov
	rm -f bin/auto_corrs
	rm -f bin/llh_scan

	rm -f build/*.o
	rm -f lib/libantinullh.a
	rm -f $(PREFIX)/make_pdfs
	rm -f $(PREFIX)/make_plots
	rm -f $(PREFIX)/make_trees
	rm -f $(PREFIX)/fit_dataset
	rm -f $(PREFIX)/build_asimov
	rm -f $(PREFIX)/auto_corrs
	rm -f $(PREFIX)/llh_scan
