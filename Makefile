CXX := g++
ROOT_FLAGS := `root-config --cflags --libs `
G4_FLAGS := `geant4-config --libs`
GSL_FLAGS := -lgsl

RAT_ROOT := $(RATROOT)
RAT_INC := $(RAT_ROOT)/include
RAT_EXTRN_INC := $(RAT_INC)/external
RAT_LIB_DIR := $(RAT_ROOT)/lib
RAT_LIB_NAME := RATEvent_Linux

OXSX_ROOT := $(OXO_DIR)
OXSX_INC := $(OXSX_ROOT)/include
OXSX_LIB_DIR := $(OXSX_ROOT)/build 
OXSX_LIB_NAME := oxsx

H5_LIBS = hdf5_hl_cpp hdf5_cpp hdf5_hl hdf5

SRC_FILES := $(wildcard src/*.cc src/*/*.cc )
OBJ_FILES := $(addprefix build/, $(notdir $(SRC_FILES:.cc=.o)))

INC_DIRS := $(wildcard src/*/ ) $(OXSX_INC) $(RAT_EXTRN_INC) $(RAT_INC)
INCLUDES := $(addprefix -I,$(INC_DIRS))

LIB_DIR := lib
LIB_NAME := antinullh
LIB := $(LIB_DIR)/lib$(LIB_NAME).a

LIB_DIRS := $(LIB_DIR) $(OXSX_LIB_DIR) $(RAT_LIB_DIR)
LIBRARYDIRS := $(addprefix -L,$(LIB_DIRS))

LIB_NAMES := $(LIB_NAME) $(OXSX_LIB_NAME) $(RAT_LIB_NAME) $(H5_LIBS)
LIBRARYNAMES := $(addprefix -l,$(LIB_NAMES))

all: bin/prune_trees bin/llh_scan bin/make_reactor_json bin/make_osc_grids bin/compare_osc_grids #bin/make_plots bin/fit_dataset bin/auto_corrs

bin/fit_dataset: exec/fit_dataset.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/fit_dataset.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) ${GSL_FLAGS} -larmadillo -lMinuit2 -o $@

bin/make_reactor_json: exec/make_reactor_json.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/make_reactor_json.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) -larmadillo -o $@

bin/make_osc_grids: exec/make_osc_grids.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/make_osc_grids.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) -o $@

bin/compare_osc_grids: exec/compare_osc_grids.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/compare_osc_grids.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) -larmadillo -o $@

bin/prune_trees: exec/prune_trees.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/prune_trees.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) -larmadillo -o $@

bin/auto_corrs: exec/auto_corrs.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/auto_corrs.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) ${GSL_FLAGS} -larmadillo -o $@

bin/llh_scan: exec/llh_scan.cc $(LIB)
	mkdir -p bin
	$(CXX)  exec/llh_scan.cc $(INCLUDES) -w $(LIBRARYDIRS) $(LIBRARYNAMES) $(ROOT_FLAGS) $(G4_FLAGS) ${GSL_FLAGS} -larmadillo -o $@

$(LIB) : $(OBJ_FILES)
	mkdir -p $(LIB_DIR)
	ar rcs  $@ $^

build/%.o : src/*/%.cc
	mkdir -p build
	$(CXX) -c -w $< $(INCLUDES) -w $(ROOT_FLAGS) $(G4_FLAGS) -o $@

clean:
	rm -f bin/make_plots
	rm -f bin/prune_trees
	rm -f bin/make_reactor_json
	rm -f bin/fit_dataset
	rm -f bin/auto_corrs
	rm -f bin/llh_scan
	rm -f bin/compare_osc_grids
	rm -f build/*.o
	rm -f lib/libantinullh.a
