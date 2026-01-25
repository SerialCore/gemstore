OBJ_DIR=obj/
SRC_DIR=src/
INCLUDE_DIR=include/

CFLAGS = -I $(INCLUDE_DIR) -lm

EXCUTEABLE=gemstore

include $(SRC_DIR)Makefile

.PHONY:

all: $(OBJ_DIR) $(OBJECT)
	gcc	$(OBJECT) -o $(EXCUTEABLE) $(CFLAGS)

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

run:
	./$(EXCUTEABLE)

clean:
	rm -rf $(OBJ_DIR)
	rm -f $(EXCUTEABLE)

install:
	sudo cp $(EXCUTEABLE) /usr/local/bin/$(EXCUTEABLE)

uninstall:
	sudo rm /usr/local/bin/$(EXCUTEABLE)