CFLAGS=-Wall -O3

all: flow

flow:flow.c
	$(CC) $(CFLAGS) flow.c -lm -o flow
