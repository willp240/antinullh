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
LIB_NAME=bbfit

LIB=$(LIB_DIR)/lib$(LIB_NAME).a

all: bin/make_pdfs bin/make_plots bin/make_trees bin/split_data bin/fit_dataset bin/up_count_lim bin/build_azimov bin/split_half bin/sum_pdfs bin/sum_pdfs_3d bin/smooth_pdfs bin/slice_pdfs bin/make_plots bin/llh_scan bin/auto_corrs

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

bin/split_data: split_data.cc $(LIB)
	mkdir -p bin
	$(CXX)  split_data.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

bin/split_half: split_half.cc $(LIB)
	mkdir -p bin
	$(CXX)  split_half.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) $(GSL_FLAGS)  -larmadillo -o $@

bin/make_trees: make_trees.cc $(LIB)
	mkdir -p bin
	$(CXX)  make_trees.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) -larmadillo -o $@

bin/build_azimov: build_azimov.cc $(LIB)
	mkdir -p bin
	$(CXX)  build_azimov.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) ${GSL_FLAGS} -larmadillo -o $@

bin/sum_pdfs: sum_pdfs.cc $(LIB)
	mkdir -p bin
	$(CXX)  sum_pdfs.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@


bin/sum_pdfs_3d: sum_pdfs_3d.cc $(LIB)
	mkdir -p bin
	$(CXX)  sum_pdfs_3d.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

bin/smooth_pdfs: smooth_pdfs.cc $(LIB)
	mkdir -p bin
	$(CXX)  smooth_pdfs.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

bin/slice_pdfs: slice_pdfs.cc $(LIB)
	mkdir -p bin
	$(CXX)  slice_pdfs.cc -I$(INC_DIR) -I$(OXSX_INC) -w -L$(LIB_DIR) -L$(OXSX_LIB_DIR) -l$(LIB_NAME) -l$(OXSX_LIB_NAME)  $(ROOT_FLAGS) $(G4_FLAGS) $(H5_LIBS) -larmadillo -o $@

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
	ln -sf `readlink -f bin/split_data` $(PREFIX)
	ln -sf `readlink -f bin/split_half` $(PREFIX)
	ln -sf `readlink -f bin/fit_dataset` $(PREFIX)
	ln -sf `readlink -f bin/fit_dataset_batch` $(PREFIX)
	ln -sf `readlink -f bin/up_count_lim` $(PREFIX)
	ln -sf `readlink -f bin/build_azimov` $(PREFIX)
	ln -sf `readlink -f bin/sum_pdfs` $(PREFIX)
	ln -sf `readlink -f bin/sum_pdfs_3d` $(PREFIX)
	ln -sf `readlink -f bin/smooth_pdfs` $(PREFIX)
	ln -sf `readlink -f bin/slice_pdfs` $(PREFIX)
	ln -sf `readlink -f bin/auto_corrs` $(PREFIX)
	ln -sf `readlink -f bin/llh_scan` $(PREFIX)
	chmod +x bin/make_pdfs
	chmod +x bin/make_plots
	chmod +x bin/make_trees
	chmod +x bin/split_data
	chmod +x bin/split_half
	chmod +x bin/fit_dataset
	chmod +x bin/fit_dataset_batch
	chmod +x bin/up_count_lim
	chmod +x bin/build_azimov
	chmod +x bin/sum_pdfs
	chmod +x bin/sum_pdfs_3d	
	chmod +x bin/smooth_pdfs
	chmod +x bin/slice_pdfs
	chmod +x bin/auto_corrs
	chmod +x bin/llh_scann

clean:
	rm -f bin/make_pdfs
	rm -f bin/make_plots
	rm -f bin/make_trees
	rm -f bin/split_data
	rm -f bin/fit_dataset
	rm -f bin/up_count_lim
	rm -f bin/build_azimov
	rm -f bin/split_half
	rm -f bin/sum_pdfs
	rm -f bin/sum_pdfs_3d
	rm -f bin/smooth_pdfs
	rm -f bin/slice_pdfs
	rm -f bin/auto_corrs
	rm -f bin/llh_scan

	rm -f build/*.o
	rm -f lib/libbbfit.a
	rm -f $(PREFIX)/make_pdfs
	rm -f $(PREFIX)/make_plots
	rm -f $(PREFIX)/make_trees
	rm -f $(PREFIX)/split_data
	rm -f $(PREFIX)/fit_dataset
	rm -f $(PREFIX)/up_count_lim
	rm -f $(PREFIX)/split_half
	rm -f $(PREFIX)/build_azimov
	rm -f $(PREFIX)/sum_pdfs
	rm -f $(PREFIX)/sum_pdfs_3d
	rm -f $(PREFIX)/smooth_pdfs
	rm -f $(PREFIX)/slice_pdfs
	rm -f $(PREFIX)/auto_corrs
	rm -f $(PREFIX)/llh_scan
