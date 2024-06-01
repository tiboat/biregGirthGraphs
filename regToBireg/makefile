compiler=gcc
flags= -lm -mcmodel=medium -O3 -march=native #-fsanitize=address

# You can adapt DMAXN according to which maximum graph order you would like to perform regToBireg.
# E.g. if the max order would be 500, then setting DMAXN to 512 would increase the performance.
regToBireg: regToBireg.c bitset
	$(compiler) -g -o regToBiregExec regToBireg.c bitset/bitset.c ../nauty/gtools.c ../nauty/nautyL.a -DMAXN=2048 -DWORDSIZE=64 $(flags)

