# Simple Makefile for CS257 coursework.
CC=gcc
LIBS=-lm -lGL -lglut
MACLIBS=-lm -framework OpenGL -framework GLUT

.PHONY: all linux apple
all: linux

linux:
	$(CC) -o cs257 -O2 -msse3 -fopenmp cs257.c $(LIBS) -std=c99

apple:
	$(CC) -o cs257 -O2 -msse3 -fopenmp cs257.c $(MACLIBS) -std=c99
