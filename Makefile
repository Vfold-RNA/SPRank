.PHONY: all sprank sprank-rf util clear

all:
	$(MAKE) sprank
	$(MAKE) sprank-rf
	$(MAKE) util
	chmod +x ./bin/*

sprank:
	mkdir -p bin
	g++ ./src/main.cpp ./src/vec3d.cpp ./src/mol.cpp -o ./bin/sprank -I ./src/ -O3 -std=c++11 -static -static-libgcc -static-libstdc++ -lstdc++fs

sprank-rf:
	cp ./random-forest/sprank-rf ./bin/sprank-rf

util:
	cp ./util/ambertools_prepare ./bin/
	cp ./util/ambertools_prepare_cpd ./bin/
	cp ./util/ambertools_prepare_rec ./bin/
	cp ./util/check_atom_order ./bin/

clear:
	rm -f bin/sprank
	rm -f bin/sprank-rf
	rm -r bin/ambertools_prepare
	rm -r bin/ambertools_prepare_cpd
	rm -r bin/ambertools_prepare_rec
	rm -r bin/check_atom_order
