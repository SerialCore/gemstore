
OBJECT += $(patsubst $(SRC_DIR)%.c, $(OBJ_DIR)%.o, $(wildcard $(SRC_DIR)*.c))

$(OBJ_DIR)%.o: $(SRC_DIR)%.c
	gcc -c $< -o $@ $(CFLAGS)