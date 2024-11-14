compiler=gcc
flags= -lm -O4 -march=native -no-pie


all: regToBireg regToBiregAddEdge

regToBireg: regToBireg.c bitset
	$(compiler) -g -o regToBiregExec regToBireg.c bitset/bitset.c ../nauty/gtools.c ../nauty/nauty.a $(flags)

regToBiregAddEdge: regToBiregAddEdge.c
	$(compiler) -g -o regToBiregAddEdgeExec regToBiregAddEdge.c ../nauty/gtools.c ../nauty/nauty.a $(flags)

clean:
	rm regToBiregExec regToBiregAddEdgeExec
