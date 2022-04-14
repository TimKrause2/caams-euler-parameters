CPPFLAGS= -ggdb `pkg-config --cflags eigen3`
CC=c++
LDLIBS= -lGL -lglut
LDFLAGS=-ggdb


n21test:n21test.o caams.o body.o constraint.o forces.o mbsystem.o

n21test.o:n21test.cpp

mb_gyro:mb_gyro.o caams.o body.o constraint.o forces.o mbsystem.o

mb_gyro.o:mb_gyro.cpp

caams.o:caams.cpp

body.o:body.cpp

eigen_test:eigen_test.cpp

constraint.o:constraint.cpp

forces.o:forces.cpp

mbsystem.o:mbsystem.cpp
