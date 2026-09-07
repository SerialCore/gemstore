CC = gcc
CPPC = g++
CFLAGS = -Wall -O2
CPPFLAGS = -Wall -O2 -DMATHCORE_STANDALONE -DROOT_Math_VecTypes

ifdef USE_LAPACKE
    CFLAGS += -DLAPACKE
    LDFLAGS += -llapacke -lopenblas
endif

INC_DIR = include/
OBJ_DIR = obj/
SRC_DIR = src/
BASIS_DIR = src/basis/
MATH_DIR = src/math/
MODEL_DIR = src/model/
PARAM_DIR = src/param/
MINUIT2_DIR = lib/Minuit2
MINUIT2_BUILD_DIR = $(MINUIT2_DIR)/build
MINUIT2_INC_DIR = $(MINUIT2_DIR)/inc
MINUIT2_LIB_DIR = $(MINUIT2_BUILD_DIR)/lib
MINUIT2_STAMP = $(MINUIT2_BUILD_DIR)/.built-stamp

EXCUTEABLE = gemstore
LIBMINUIT2 = $(MINUIT2_LIB_DIR)/libMinuit2.a
LIBMINUIT2MATH = $(MINUIT2_LIB_DIR)/libMinuit2Math.a
LDFLAGS = $(LIBMINUIT2) $(LIBMINUIT2MATH) -lm

INCLUDES = -I$(INC_DIR)
INCLUDES_CPP = $(INCLUDES) -I$(MINUIT2_INC_DIR)

C_SOURCES = $(wildcard $(SRC_DIR)*.c $(BASIS_DIR)*.c $(MATH_DIR)*.c $(MODEL_DIR)*.c $(PARAM_DIR)*.c)
CPP_SOURCES = $(wildcard $(PARAM_DIR)*.cc)

OBJECTS = $(patsubst %.c, $(OBJ_DIR)%.o, $(C_SOURCES)) \
          $(patsubst %.cc, $(OBJ_DIR)%.o, $(CPP_SOURCES))

all: $(LIBMINUIT2) $(LIBMINUIT2MATH) $(EXCUTEABLE)

$(EXCUTEABLE): $(OBJECTS) $(LIBMINUIT2) $(LIBMINUIT2MATH)
	${CPPC} $^ -o $@ $(LDFLAGS)


$(LIBMINUIT2) $(LIBMINUIT2MATH): $(MINUIT2_STAMP)

$(MINUIT2_STAMP):
	@test -f "$(MINUIT2_DIR)/CMakeLists.txt" || (echo "Minuit2 submodule is not initialized. Run: git submodule update --init --recursive" && false)
	cmake -S "$(MINUIT2_DIR)" -B "$(MINUIT2_BUILD_DIR)" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
	cmake --build "$(MINUIT2_BUILD_DIR)" --target Minuit2 Minuit2Math
	@touch "$@"

$(OBJ_DIR)%.o: %.c
	@mkdir -p $(@D)
	${CC} $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR)%.o: %.cc
	@mkdir -p $(@D)
	${CPPC} $(CPPFLAGS) $(INCLUDES_CPP) -c $< -o $@

run: all
	./$(EXCUTEABLE)

clean:
	rm -rf $(OBJ_DIR) $(EXCUTEABLE)
	rm -rf $(MINUIT2_BUILD_DIR)

install: all
	cp $(EXCUTEABLE) ~/.local/bin/

uninstall:
	rm -f ~/.local/bin/$(EXCUTEABLE)

print-objects:
	@echo "OBJECTS = $(OBJECTS)"
	@echo "LIBMINUIT2 = $(LIBMINUIT2)"
	@echo "LIBMINUIT2MATH = $(LIBMINUIT2MATH)"

.PHONY: all run clean install uninstall print-objects
