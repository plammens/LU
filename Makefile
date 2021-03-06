############################# User-defined constants ################################

##### Directory names #####

INCLUDE_DIR := include
SRC_DIR := src
LIB_ROOT_DIR := libs
TEST_DIR := test

# directory for build files (binaries, object files, etc.):
BUILD_DIR := build


##### Configuration #####

# name of main executable:
MAIN_NAME := main

# test suite and test executable name:
TEST_SUITE := doctest
TEST_NAME := test




############################# Macros ################################


##### Directory names #####

TEST_TEXT_DIR = $(TEST_DIR)/text
BIN_DIR = $(BUILD_DIR)/bin
OBJ_DIR = $(BUILD_DIR)/obj

DEP_DIR = $(BUILD_DIR)/depend
LIB_INCLUDE_DIR = $(LIB_ROOT_DIR)/include
LIB_FILE_DIR = $(LIB_ROOT_DIR)/lib


##### file paths #####

MAIN_EXE = $(BIN_DIR)/$(MAIN_NAME).x
MAIN_OBJ = $(OBJ_DIR)/$(MAIN_NAME).o

TEST_EXE = $(BIN_DIR)/$(TEST_NAME).x



##### Compiler options and flags ######

CXX = g++
CXXFLAGS = -std=c++11 -O2

CXX_COMPILE_FLAGS = $(CXXFLAGS) -I $(INCLUDE_DIR) -I $(LIB_INCLUDE_DIR)
CXX_LINK_FLAGS = $(CXXFLAGS)

CXX_TEST_COMPILE_FLAGS = $(CXX_COMPILE_FLAGS) -I $(TEST_DIR)/$(INCLUDE_DIR)
CXX_TEST_LINK_FLAGS = $(CXX_LINK_FLAGS)


##### Auto-detected files and paths #####

# Directory search paths:
vpath %.h $(shell find $(INCLUDE_DIR) -type d)
vpath %.cpp $(shell find $(SRC_DIR) -type d)
vpath test%.cpp $(shell find $(TEST_DIR)/$(SRC_DIR) -type d)

# sources and objects:
sources := $(shell find $(SRC_DIR) -type f -name '*.cpp' ! -name '$(MAIN_NAME)*')  # sources excluding main
objects := $(patsubst %.cpp,$(OBJ_DIR)/%.o, $(notdir $(sources)))  # object files

# test sources (prefixed with 'test'):
test_sources := $(shell find $(TEST_DIR)/$(SRC_DIR) -type f -name 'test*.cpp')
test_objects := $(patsubst %.cpp,$(OBJ_DIR)/%.o, $(notdir $(test_sources)))

# dependency files for automatic Makefile rule prerequisites
depends := $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$(objects) $(test_objects))



################################## Rules ######################################

############### Phony rules ###############

.PHONY: all build build-test debug libs run test \
		clean clean-build clean-out docs view-docs examples \
		.pre-build .pre-lib .pre-build-test


all: build build-test examples


build: .pre-build $(MAIN_EXE)
	@printf "\e[1mDone building main.\e[0m\n\n"

build-test: debug .pre-build-test $(TEST_EXE)
	@printf "\e[1mDone building tests.\e[0m\n\n"

debug: CXXFLAGS += -Wall -Og
debug: build


run: build
	@echo
	@printf "\e[1mExecuting main program...\e[0m\n"
	$(MAIN_EXE) $(ARGS)
	@echo
	@echo

test: build-test
	@echo
	@printf "\e[1mStarting $(TEST_SUITE)...\e[0m\n"
	$(TEST_EXE) $(ARGS)
	@echo
	@echo


clean: clean-build clean-out

clean-build:
	rm -rf ./$(BUILD_DIR)

clean-out:
	rm -rf SOLUTION*.DAT


# Just some text:

.pre-build:
	@printf "\e[1mBuilding main program...\e[0m\n"

.pre-build-test:
	@printf "\e[1mBuilding tests...\e[0m\n"



########## Main project rules #############

-include $(depends)  # include dependecy relationships

# This rule links all object files together and outputs the main executable.
$(MAIN_EXE): $(objects) $(MAIN_OBJ) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(CXX_LINK_FLAGS)

# This rule compiles source files into their corresponding object file.
# As a side effect of compilation we generate a dependency file readable
# by make (with the `-MMD -MF` flags)
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR) $(DEP_DIR)
	$(CXX) -c $< -o $@ $(CXX_COMPILE_FLAGS) -MMD -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@)


# Create build directories if nonexistent:

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(DEP_DIR):
	mkdir -p $(DEP_DIR)



############### Test rules ################

# Compile and link all test objects (together with the project's objects), outputting an executable:
$(TEST_EXE): $(objects) $(test_objects) | $(BIN_DIR)
	$(CXX) $^ -o $@ $(CXX_TEST_LINK_FLAGS)

# This rule compiles test source files into their corresponding object file.
# As a side effect of compilation we generate a dependency file.
$(OBJ_DIR)/test%.o: test%.cpp | $(OBJ_DIR) $(DEP_DIR)
	$(CXX) -c $< -o $@ $(CXX_TEST_COMPILE_FLAGS) -MMD -MF $(patsubst $(OBJ_DIR)/%.o,$(DEP_DIR)/%.d,$@)

