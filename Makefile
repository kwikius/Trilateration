
INCLUDES = -I/home/andy/cpp/projects/quan-trunk/

CXXFLAGS = -fconcepts -std=c++17 -O2

objects = trilateration_transform_matrix_minimal.o

all : test.exe

CXX = g++-7

test.exe : ${objects}
	$(CXX) -o $@  $<

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@
	$(CXX) $(CXXFLAGS) $(INCLUDES) -S $< -o main.asm

clean :
	rm -f *.o *.exe *.asm

