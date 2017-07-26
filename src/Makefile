CC=gcc
CXX=g++
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: CHROMEISTER



CHROMEISTER: CHROMEISTER.c
	$(CC) $(CFLAGS) alignmentFunctions.c -lpthread -lm commonFunctions.c -lm CHROMEISTER.c -lpthread -lm -o $(BIN)/CHROMEISTER

clean:
	rm -rf $(BIN)/CHROMEISTER
