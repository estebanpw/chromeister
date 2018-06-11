CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: lowcompfilter



lowcompfilter: lowcompfilter.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lm commonFunctions.c -lm lowcompfilter.c -lm -o $(BIN)/lowcompfilter

clean:
	rm -rf $(BIN)/lowcompfilter
