CPPC = g++
CPPFLAGS =
INCLUDES = -I $(INC_DIR)

OBJ_DIR = obj/
SRC_DIR = src/
INC_DIR = include/

TARGET = libMinuit2.a

OBJECTS = $(patsubst %.cxx, $(OBJ_DIR)%.o, $(wildcard $(SRC_DIR)*.cxx))

all: $(TARGET)

$(TARGET): $(OBJECTS)
	ar rcs $@ $^

$(OBJ_DIR)%.o: %.cxx
	@mkdir -p $(@D)
	$(CPPC) $(CPPFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(TARGET)

print-objects:
	@echo "OBJECTS = $(OBJECTS)"

.PHONY: all clean print-objects