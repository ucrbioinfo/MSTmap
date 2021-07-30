CC = g++

all: mstmap

mstmap: main.o linkage_group_DH.o genetic_map_DH.o  linkage_group_RIL.o genetic_map_RIL.o MSTOpt.o
	$(CC) -O3 -o mstmap main.o linkage_group_DH.o genetic_map_DH.o linkage_group_RIL.o genetic_map_RIL.o MSTOpt.o

main.o: main.cpp constants.h
	$(CC) -O3 -o main.o -c -g main.cpp

MSTOpt.o: MSTOpt.cpp MSTOpt.h
	$(CC) -O3 -o MSTOpt.o -c -g MSTOpt.cpp

genetic_map_DH.o: genetic_map_DH.cpp genetic_map_DH.h
	$(CC) -O3 -o genetic_map_DH.o -c -g genetic_map_DH.cpp 

linkage_group_DH.o : linkage_group_DH.h linkage_group_DH.cpp
	$(CC) -O3 -o linkage_group_DH.o  -c -g linkage_group_DH.cpp

genetic_map_RIL.o: genetic_map_RIL.cpp genetic_map_RIL.h
	$(CC) -O3 -o genetic_map_RIL.o -c -g genetic_map_RIL.cpp 

linkage_group_RIL.o : linkage_group_RIL.h linkage_group_RIL.cpp
	$(CC) -O3 -o linkage_group_RIL.o  -c -g linkage_group_RIL.cpp

clean: 
	rm -rf *.o; rm -rf mstmap; rm -rf *~; rm -rf *.out
