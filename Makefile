all:
	gcc src/main.c src/cg.c src/fileio.c -I include/ -lm -o gemstore

clean:
	rm gemstore
