all:
	gcc -o gemstore src/basis.c src/debug.c src/eigen.c src/fileio.c src/inteCenV.c src/main.c src/matrix.c src/mfi.c src/parallel.c src/spin.c src/sumckdk.c src/vtype.c -I include/ -lm

clean:
	rm gemstore