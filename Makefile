CPPFLAGS=-O2
CC=g++
LDLIBS=-lGL -lglut

all:3_piece_mobile test 3_piece_mobile_integrate 3_piece_mobile_gyro 3_piece_mobile_opengl unconstrained

3_piece_mobile_opengl:3_piece_mobile_opengl.o matrix.o caams.o

3_piece_mobile_gyro:3_piece_mobile_gyro.o matrix.o caams.o

3_piece_mobile_integrate:3_piece_mobile_integrate.o matrix.o caams.o

3_piece_mobile:3_piece_mobile.o matrix.o caams.o

test:test.o matrix.o caams.o

unconstrained:unconstrained.o matrix.o caams.o

unconstrained.o:unconstrained.cpp

test.o:test.cpp

matrix.o:matrix.cpp

caams.o:caams.cpp

3_piece_mobile.o:3_piece_mobile.cpp

3_piece_mobile_integrate.o:3_piece_mobile_integrate.cpp

3_piece_mobile_gyro.o:3_piece_mobile_gyro.cpp

3_piece_mobile_opengl.o:3_piece_mobile_opengl.cpp
