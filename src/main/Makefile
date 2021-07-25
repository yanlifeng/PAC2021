include ../../Makefile.inc
BIN_MODULE := $(patsubst %.c, %_c, $(wildcard *.c))
BIN_MODULE += $(patsubst %.cc, %_cc, $(wildcard *.cc))
BIN_MODULE += $(patsubst %.cpp, %_cpp, $(wildcard *.cpp))
BIN_MODULE += $(patsubst %.cxx, %_cxx, $(wildcard *.cxx))

define link-cmd
	@rm -rf $(BIN_MODULE)
	$(CC) $(CFLAGS) -I$(INC_PATH) -I$(EXTERNAL_INC_PATH) -I$(EXTERNAL_INC_EIGEN_PATH) -L$(LIB_PATH) -L$(EXTERNAL_LIB_PATH) -o $(addprefix ../../bin/, $@) $< $(LIBS)
endef
#This command is used to generate -g version
define link-cmd-g
	@rm -rf $(BIN_MODULE)
	$(CC) $(CFLAGS_g) -I$(INC_PATH) -I$(EXTERNAL_INC_PATH) -I$(EXTERNAL_INC_EIGEN_PATH) -L$(LIB_PATH) -L$(EXTERNAL_LIB_PATH) -o $(addprefix ../../bin/, $@)_g $< $(LIBS)
endef
VPATH = ../../bin
all:$(BIN_MODULE)
%_c:%.c
	$(link-cmd)
%_cc:%.cc
	$(link-cmd)
%_cpp:%.cpp
	$(link-cmd)
	$(link-cmd-g)
%_cxx:%.cxx
	$(link-cmd)
.phony:clean
clean:
	rm -rf $(BIN_MODULE)

