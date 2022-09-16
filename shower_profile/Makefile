all: main

main: main.o shower_creator.o ice_profile.o shower_1D.o charge_excess_profile.o
	g++ -o main main.o shower_creator.o ice_profile.o shower_1D.o charge_excess_profile.o

	
main.o: main.cpp
	g++ -c main.cpp


shower_1D.o: shower_1D.cpp shower_1D.h
	g++ -c shower_1D.cpp

ice_profile.o: ice_profile.cpp ice_profile.h
	g++ -c ice_profile.cpp

shower_creator.o: shower_creator.cpp shower_creator.h ice_profile.h shower_1D.h
	g++ -c shower_creator.cpp

charge_excess_profile.o: charge_excess_profile.cpp
	g++ -c charge_excess_profile.cpp

clean:
	rm *.o main