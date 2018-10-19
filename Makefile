CC=gcc
LIBS=-lm
THREAD_LIBS=-lpthread -lrt
PRE_CFLAGS=-pg -O0
CFLAGS=$(PRE_CFLAGS) -std=gnu99 -Wall

default: nbody

nbody: src/nbody.c thread_control.o thread_barrier.o
	$(CC) -o $@ $^ $(LIBS) $(CFLAGS) $(THREAD_LIBS)

thread_barrier.o: src/thread_barrier.c src/thread_barrier.h
	$(CC) -c -o $@ $< $(LIBS) $(CFLAGS)

thread_control.o: src/thread_control.c src/thread_control.h thread_barrier.o
	$(CC) -c -o $@ $< $(LIBS) $(CFLAGS)

