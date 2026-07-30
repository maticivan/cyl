CXX = $(shell ls /opt/homebrew/bin/g++-* | sort -V | tail -1)
 
all:
	g++ -I /usr/local/include -c cylinders.cpp -std=c++23 -lstdc++fs -pthread -fmax-errors=1
	g++ -o cyl cylinders.o -std=c++23 -lstdc++fs -pthread

forMac:
	$(CXX) -I /usr/local/include -c cylinders.cpp -std=c++23 -pthread -DIOF_NO_SPAWN -fmax-errors=1
	$(CXX) -o cyl cylinders.o -std=c++23 -lstdc++fs -pthread
