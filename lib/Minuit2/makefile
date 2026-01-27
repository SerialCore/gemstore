CCPP = g++
CPP_ARGS =

OBJ_DIR = obj
MINUIT2_SRC_DIR = src
MINUIT2_INC_DIR = include

FILE_CXX := $(wildcard ./$(MINUIT2_SRC_DIR)/*.cxx)
OBJS_CXX  = $(addprefix ./$(OBJ_DIR)/,$(notdir $(FILE_CXX:.cxx=.o))) 

all: libMinuit2.a a.out

$(OBJ_DIR)/%.o: $(MINUIT2_SRC_DIR)/%.cxx
	mkdir -p $(OBJ_DIR) && $(CCPP)$(CCPP_ARGS) -c $< -o $@ -I $(MINUIT2_INC_DIR)

libMinuit2.a: $(OBJS_CXX)
	ar rcs libMinuit2.a $(OBJS_CXX)

a.out: find_point.cpp libMinuit2.a
	$(CCPP) $(CCPP_ARGS) -o a.out find_point.cpp libMinuit2.a -I $(MINUIT2_INC_DIR)

clean:
	rm -rf libMinuit2.a a.out $(OBJ_DIR)

print-x:
	@echo $(OBJS_CXX)
