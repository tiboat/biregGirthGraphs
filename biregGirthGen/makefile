compiler=gcc
flags= -lm #-O3 -march=native #-fsanitize=address

## Compile commands for up to 64 vertices
biregGirth: biregGirthGen.c bitset
	$(compiler) -DUSE_64_BIT -g -o biregGirthGen biregGirthGen.c ../nauty/nautyL1.a -DMAXN=WORDSIZE $(flags)

biregGirth-profiler: biregGirthGen.c bitset
	$(compiler) -DUSE_64_BIT -o biregGirthGen-profiler biregGirthGen.c ../nauty/nautyL1.a -DMAXN=WORDSIZE -pg $(flags)


## Compile commands for up to 128 vertices

# There are two different implementations of the 128-bit version. The array version generally performs faster. (source: GenK2)
biregGirth-128: biregGirthGen.c bitset
	$(compiler) -DUSE_128_BIT -g -o biregGirthGen-128 biregGirthGen.c ../nauty/nautyL.a -DMAXN=128 $(flags)

biregGirth-128-profiler: biregGirthGen.c bitset
	$(compiler) -DUSE_128_BIT -o biregGirthGen-128-profiler biregGirthGen.c ../nauty/nautyL.a -DMAXN=128 -pg $(flags)

biregGirth-128-array: biregGirthGen.c bitset
	$(compiler) -DUSE_128_BIT_ARRAY -g -o biregGirthGen-128-array biregGirthGen.c ../nauty/nautyL.a -DMAXN=128 $(flags)

biregGirth-128-array-profiler: biregGirthGen.c bitset
	$(compiler) -DUSE_128_BIT_ARRAY -o biregGirthGen-128-array-profiler biregGirthGen.c ../nauty/nautyL.a -DMAXN=128 -pg $(flags)



clean:
	rm -f biregGirth biregGirth-profiler biregGirth-128 biregGirth-128-profiler biregGirth-128-array biregGirth-128-array-profiler

