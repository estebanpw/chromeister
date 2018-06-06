CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: CHROMEISTER dipeptides



CHROMEISTER: CHROMEISTER.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm CHROMEISTER.c -lm -o $(BIN)/CHROMEISTER

dipeptides: dipeptides.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm dipeptides.c -lm -o $(BIN)/dipeptides

clean:
	rm -rf $(BIN)/CHROMEISTER $(BIN)/dipeptides
