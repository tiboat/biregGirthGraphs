# calcMaxAtMinDist


This directory contains code to calculate the maximum (or a large) amount of vertices or edges at minimum pairwise distance $\left\lceil g/2\right\rceil$ in an $(r,g)$-graph ($=$ $r$-regular graph of girth $g$).

## Data

### Lists of regular graphs

The file `regularGraphs.s6` contains $(r,g)$-graphs that we found online (or constructed ourselves). The following websites contain the $(r,g)$-graphs we used.

 - http://school.maths.uwa.edu.au/~gordon/remote/cages/index.html (access via the Wayback machine)
 - http://ginger.indstate.edu/ge/CAGES (access via the Wayback machine)
 - http://cs.indstate.edu/ge/SmallCages/index.html
 - https://aeb.win.tue.nl/graphs/cages/cages.html
 - https://github.com/FransdR/Record-graphs-girth-7-and-8

We programmed the construction for the $(5,11)$-graph of order 2688, as described in [1]. The Python script `oddGirthConstr.py` constructs and outputs this graph in sparse6 format. The construction starts from a $(5,12)$ graph, which is seperately available in `data/5_12.s6`.


### Best number of vertices/edges found

The file `distRegularGraphs.txt` contains a summary of the largest number of vertices and edges we found at minimum pairwise distance $\left\lceil g/2\right\rceil$. These numbers are not all (necessarily) maximal, since for some of these, approximate calculations where used (as explained underneath).


## Code

### Compilation


First configure and compile nauty. Go to the nauty directory (`biregGirthGraphs/nauty/`) and execute the following two commands.

```
./configure
make
```

Then, since this code uses OpenMP, you might want to install this first (if you haven't isntalled this yet).

Then compile `calcMaxAtMinDist.c` by using the `make` command, which creates the executable `calcMaxAtMinDistExec`.
```
make
```

### Execution

To execute `calcMaxAtMinDistExec`, feed it with a list of $r$-regular graphs via stdin and give the executable the proper arguments.

```
cat <list_file> | ./calcMaxAtMinDistExec [-v|-e|-g]
````

Here `<list_file>` is a list of $r$-regular graphs in sparse6 or graph6 format and `-v`, `-e` and `-g` are optional flages with the following meaning.

 * `-v`: output the most amount of vertices at the minimum pairwise distance found so far,
 * `-e`: output the most amount of edges at the minimum pairwise distance found so far, and
 * `-g`: applies a greedy search for finding a number of vertices and edges at the pairwise minimum distance. This number is not necessarily optimal.

When no argument is given, it computes the maximum amount of vertices and edges exactly and only shows the output when the computation is finished.

These optional arguments are useful for large graphs where finding the exact amount of vertices and edges at the minimum pairwise distance can be too expensive.


#### Example no arguments
```
cat data/regularGraphs.s6 | ./calcMaxAtMinDistExec
````
This goes through the list of regular graphs in the file `regularGraphs.s6` and prints the maximum amount of vertices and edges at the minimum pairwise distance $\left\lceil g/2\right\rceil$. The output looks like this.

```
n=10 (3,5) diam=2
num vertices at pairwise distance at least 3: 1
num edges at pairwise distance at least 3: 1
n=14 (3,6) diam=3
num vertices at pairwise distance at least 3: 2
num edges at pairwise distance at least 3: 1
n=19 (4,5) diam=3
num vertices at pairwise distance at least 3: 3
num edges at pairwise distance at least 3: 2
...
```



#### Example `-g`

```
cat data/regularGraphs.s6 | ./calcMaxAtMinDistExec -g
````
This goes through the list of regular graphs in the file `regularGraphs.s6` and prints an amount of vertices and edges at the minimum pairwise distance $\left\lceil g/2\right\rceil$, found greedily. The output is analogous to execution with no arguments, but the numbers found might not be maximal.

```
Greedily (not exact) searching for vertices and edges at min dist
n=10 (3,5) diam=2
num vertices at pairwise distance at least 3: 1
num edges at pairwise distance at least 3: 1
n=14 (3,6) diam=3
num vertices at pairwise distance at least 3: 2
num edges at pairwise distance at least 3: 1
n=19 (4,5) diam=3
num vertices at pairwise distance at least 3: 3
num edges at pairwise distance at least 3: 2
...
```
For these three examples shown, the greedy search ends up finding the same amount of vertices and edges. Yet, this is not always the case.

#### Example `-v`
```
cat data/regularGraphs.s6 | ./calcMaxAtMinDistExec -v
````
This goes through the list of regular graphs in the file `regularGraphs.s6` and prints the largest amount of vertices at the minimum pairwise distance $\left\lceil g/2\right\rceil$ found so far. The output looks like this.

```
Printing the maximum amount of vertices at min dist found so far
n=10 (3,5) diam=2
Best maxVerticesSoFar = 1
n=14 (3,6) diam=3
Best maxVerticesSoFar = 1
Best maxVerticesSoFar = 2
n=19 (4,5) diam=3
Best maxVerticesSoFar = 1
Best maxVerticesSoFar = 2
Best maxVerticesSoFar = 3
...
```

#### Example `-e`
```
cat data/regularGraphs.s6 | ./calcMaxAtMinDistExec -e
````
This goes through the list of regular graphs in the file `regularGraphs.s6` and prints the largest amount of edges at the minimum pairwise distance $\left\lceil g/2\right\rceil$ found so far. The output looks like this.

```
Printing the maximum amount of edges at min dist found so far
n=10 (3,5) diam=2
Best maxEdgesSoFar = 1
n=14 (3,6) diam=3
Best maxEdgesSoFar = 1
n=19 (4,5) diam=3
Best maxEdgesSoFar = 1
Best maxEdgesSoFar = 2
...
```




## References

[1] G. Araujo-Pardo, On upper bounds of odd girth cages, Discrete Math. 310 (2010) 1622-1626.

