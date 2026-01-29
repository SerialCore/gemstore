CC=gcc

CFLAGS = -I $(INCLUDE_DIR) -lm

INCLUDE_DIR=include/
LIB_DIR=lib/
OBJ_DIR=obj/
SRC_DIR=src/
MODEL_DIR=src/model/
BASIS_DIR=src/basis/
NUMERICAL_DIR=src/numerical

EXCUTEABLE=gemstore

#include $(SRC_DIR)Makefile.in
#include $(MODEL_DIR)Makefile.in
#include $(BASIS_DIR)Makefile.in

OBJECT += $(patsubst $(BASIS_DIR)%.c, $(OBJ_DIR)%.o, $(wildcard $(BASIS_DIR)*.c))
OBJECT += $(patsubst $(MODEL_DIR)%.c, $(OBJ_DIR)%.o, $(wildcard $(MODEL_DIR)*.c))
OBJECT += $(patsubst $(NUMERICAL_DIR)%.c, $(OBJ_DIR)%.o, $(wildcard $(NUMERICAL_DIR)*.c))
OBJECT += $(patsubst $(SRC_DIR)%.c, $(OBJ_DIR)%.o, $(wildcard $(SRC_DIR)*.c))

all: $(OBJ_DIR) $(OBJECT)
	${CC} $(OBJECT) -o $(EXCUTEABLE) $(CFLAGS)
.PHONY: all

$(OBJ_DIR)%.o: $(BASIS_DIR)%.c
	${CC} -c $< -o $@ $(CFLAGS)

$(OBJ_DIR)%.o: $(MODEL_DIR)%.c
	${CC} -c $< -o $@ $(CFLAGS)

$(OBJ_DIR)%.o: $(NUMERICAL_DIR)%.c
	${CC} -c $< -o $@ $(CFLAGS)

$(OBJ_DIR)%.o: $(SRC_DIR)%.c
	${CC} -c $< -o $@ $(CFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

run:
	./$(EXCUTEABLE)
.PHONY: run

clean:
	rm -rf $(OBJ_DIR)
	rm -f $(EXCUTEABLE)
.PHONY: clean

install:
	sudo cp $(EXCUTEABLE) /usr/local/bin/$(EXCUTEABLE)
.PHONY: install

uninstall:
	sudo rm /usr/local/bin/$(EXCUTEABLE)
.PHONY: uninstall

print-objects:
	@echo $(OBJECT)
.PHONY: print-objects
