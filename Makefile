CC = gcc
CPPC = g++
CFLAGS = -Wall -O2
CPPFLAGS = -Wall -O2
LDFLAGS = -L$(MINUIT2_DIR) -lMinuit2 -lm

ifdef USE_LAPACKE
    CFLAGS += -DLAPACKE
    LDFLAGS += -llapacke -lopenblas
endif

INC_DIR = include/
OBJ_DIR = obj/
SRC_DIR = src/
MODEL_DIR = src/model/
BASIS_DIR = src/basis/
NUMERICAL_DIR = src/numerical/
MINUIT2_DIR = lib/Minuit2/

EXCUTEABLE = gemstore
LIBMINUIT2 = $(MINUIT2_DIR)libMinuit2.a

INCLUDES = -I$(INC_DIR)
INCLUDES_CPP = $(INCLUDES) -I$(MINUIT2_DIR)include/

C_SOURCES = $(wildcard $(SRC_DIR)*.c $(MODEL_DIR)*.c $(BASIS_DIR)*.c $(NUMERICAL_DIR)*.c)
CPP_SOURCES = $(wildcard $(SRC_DIR)*.cc)

OBJECTS = $(patsubst %.c, $(OBJ_DIR)%.o, $(C_SOURCES)) \
          $(patsubst %.cc, $(OBJ_DIR)%.o, $(CPP_SOURCES))

all: $(LIBMINUIT2) $(EXCUTEABLE)

$(EXCUTEABLE): $(OBJECTS) $(LIBMINUIT2)
	${CPPC} $^ -o $@ $(LDFLAGS)

$(LIBMINUIT2):
	$(MAKE) -C $(MINUIT2_DIR)

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
	$(MAKE) -C $(MINUIT2_DIR) clean

install: all
	sudo cp $(EXCUTEABLE) /usr/local/bin/

uninstall:
	sudo rm -f /usr/local/bin/$(EXCUTEABLE)

print-objects:
	@echo "OBJECTS = $(OBJECTS)"
	@echo "LIBMINUIT2 = $(LIBMINUIT2)"

.PHONY: all run clean install uninstall print-objects