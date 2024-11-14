compiler=gcc
flags= -lm  -march=native -no-pie -O4 -fopenmp #-fsanitize=address

calcMaxAtMinDist: calcMaxAtMinDist.c
	$(compiler) -g -o calcMaxAtMinDistExec calcMaxAtMinDist.c ../nauty/gtools.c ../nauty/nauty.a $(flags)
