sprank:
	mkdir -p bin
	g++ ./src/main.cpp ./src/vec3d.cpp ./src/mol.cpp -o ./bin/sprank -I ./src/ -O3 -std=c++11 -static -static-libgcc -static-libstdc++ -lstdc++fs
.PHONY : sprank

clear:
	rm -f bin/sprank
.PHONY : clear
