CPPFLAGS = -ggdb -O2
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

system:system.o matrix.o caams.o

system.o:system.cpp

spring:spring.o matrix.o caams.o

spring.o:spring.cpp

body.o:body.cpp

constraint.o:constraint.cpp

mbsystem.o:mbsystem.cpp

mbtest:mbtest.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

mbtest.o:mbtest.cpp

forces.o:forces.cpp

sp_sp_test:sp_sp_test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

sp_sp_test.o:sp_sp_test.cpp

rev_test:rev_test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

rev_test.o:rev_test.cpp

revp_test:revp_test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

revp_test.o:revp_test.cpp

spring_test:spring_test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

spring_test.o:spring_test.cpp

p22test:p22test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

p22test.o:p22test.cpp

n21test:n21test.o matrix.o caams.o mbsystem.o body.o constraint.o forces.o

n21test.o:n21test.cpp
