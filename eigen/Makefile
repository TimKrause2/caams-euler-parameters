CPPFLAGS= -O2 `pkg-config --cflags eigen3`
CC=c++
LDLIBS= -lGL -lglut

all:sp_sp_test n21test mb_gyro revp_test p22test chain

chain:chain.o caams.o body.o constraint.o forces.o mbsystem.o

chain.o:chain.cpp

p22test:p22test.o caams.o body.o constraint.o forces.o mbsystem.o

p22test.o:p22test.cpp

revp_test:revp_test.o caams.o body.o constraint.o forces.o mbsystem.o

revp_test.o:revp_test.cpp

sp_sp_test:sp_sp_test.o caams.o body.o constraint.o forces.o mbsystem.o

sp_sp_test.o:sp_sp_test.cpp

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
