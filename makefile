# make file for Web Data Compression and Search Assignment 1 by Felix Yuen Dao Phu

CC = g++
CFLAGS  = -g -Wall

all: program1 program2

program1: ahencode.cpp
	$(CC) $(CFLAGS) -o ahencode ahencode.cpp
	
program2: ahdecode.cpp
	$(CC) $(CFLAGS) -o ahdecode ahdecode.cpp