CPPFLAGS= -O2 `pkg-config --cflags eigen3` `sdl2-config --cflags`
CC=c++
LDLIBS= -lGL -lglut `sdl2-config --libs`

all:sp_sp_test n21test mb_gyro revp_test p22test chain vehicle

vehicle:vehicle.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

vehicle.o:vehicle.cpp

chain:chain.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

chain.o:chain.cpp

p22test:p22test.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

p22test.o:p22test.cpp

revp_test:revp_test.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

revp_test.o:revp_test.cpp

sp_sp_test:sp_sp_test.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

sp_sp_test.o:sp_sp_test.cpp

n21test:n21test.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

n21test.o:n21test.cpp

mb_gyro:mb_gyro.o caams.o body.o constraint.o forces.o mbsystem.o primitives.o

mb_gyro.o:mb_gyro.cpp

caams.o:caams.cpp

body.o:body.cpp

eigen_test:eigen_test.cpp

constraint.o:constraint.cpp

forces.o:forces.cpp

mbsystem.o:mbsystem.cpp

primitives.o:primitives.cpp

