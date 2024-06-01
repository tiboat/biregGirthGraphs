# regToBireg

This directory contains code to try to construct $(\{r,m\};g)$-graphs based on $r$-regular graphs. Given a list of $r$-regular graphs, $m$, $g$ and the range $[minN, maxN]$ in which the order of the $(\{r,m\};g)$-graph should be, the code loops over the list. For each $r$-regular graph $G$ it first checks three necessary conditions on the order, girth and diameter of $G$, to make it feasible to obtain an $(\{r,m\};g)$-graph with order $n \in [minN, maxN]$. Only if these three conditions are met, then the construction is tried. The construction deletes a limited set of edges and connects the vertices which lost an edge to either one or two added vertices, which will have degree $m$. If $m$ is even, one vertex of degree $m$ is added. If $m$ is odd, two vertices of degree $m$ are added.

## Data

### Lists of regular graphs

This directory contains compressed versions of lists or regular graphs. More precisely, 
 -  Files starting with `cubicvt` contain $3$-regular vertex-transitive graphs, where the two following numbers indicate the range of order of the upcoming files. These where converted from sparse6 to graph6 format from the repository [ cubic-vertextransitive-graphs](https://github.com/kguo-sagecode/cubic-vertextransitive-graphs). Yet, the list originates from [1,2].
 -  `cubic_symmetric` contains lists of $3$-regular symmetric graphs [3,4].
 -  `quadric` contains lists of $4$-regular graphs which are edge-transitive, arc-transitive or 2-arc-transitive [5,6,7,8,9].
 -  `penta_arc_trans` contains a list of $5$-regular arc-transitive graphs [10].
   
Unzip these files in order to use the lists inside of them.


### Obtained record graphs

The directory [found_graphs](found_graphs/) contains $(\{r,m\};g)$-graphs encoded in graph6 format of smaller order than we found in the literature. Each graph is stored in a different file formatted as `n<n>_r<r>_m<m>_g<g>.g6`, e.g. `n545_r3_m4_g14.g6`.


## Code

### Compilation

In order to execute this code, one should first compile by using the `make` command, which creates the executable `regToBiregExec`.
```
make
```

### Execution

Then one should feed the executable with a list of $r$-regular graphs via stdin and give the executable the proper arguments.
```
cat <list_file> | ./regToBiregExec <r> <m> <g> <minN> <maxN>
````
where `<list_file>` is the list of $r$-regular graphs in graph6 format and `<r>`, `<m>`, `<g>`, `<minN>` and `<maxN>` have the same meaning as described above.

#### Example
```
cat graph_lists/cubicvt502-600_graph6.g6 | ./regToBiregExec 3 4 14 328 619
````
This goes through the list of $3$-regular vertex-transitive graphs of order in $[502,600]$ in order to try to construct a $(\{3,4\};14)$-graph of order $n \in [328, 619]$. (Note that the `minN` and `maxN` values are larger than needed in this case, since the construction only adds one vertex.)

### Output
The code outputs the order, girth and diameter of the regular graphs that satisfy the three conditions. If an $(\{r,m\};g)$-graph of order $n \in [minN, maxN]$ is found by the construction, this graph is printed in graph6 format and the code does not continue trying to construct more graphs.

#### Example
Executing
```
cat graph_lists/quadric-arc-trans.g6 | ./regToBiregExec 4 7 9 281 742
```
gives the output
```
diameter=7
initial girth=9
initial n=320
diameter=6
initial girth=9
initial n=320
diameter=8
initial girth=9
initial n=360
diameter=8
initial girth=10
initial n=420
n, new girth: 422, 9
~?Ee@?G?C?? [...] ???_
```
where the complete graph6 string is trimmed to `~?Ee@?G?C?? [...] ???_`. 



## References

[1]     P. Potočnik, P. Spiga, G. Verret, Cubic vertex-transitive graphs on up to 1280 vertices, Journal of Symbolic Computation 50 (2013), 465-477. 

[2] P. Potočnik, P. Spiga, G. Verret, Bounding the order of the vertex-stabiliser in 3-valent vertex-transitive and 4-valent arc-transitive graphs, arXiv:1010.2546v1 [math.CO]. 

[3] M. Conder. Trivalent symmetric graphs on up to 768 vertices. Journal of
Combinatorial Mathematics and Combinatorial Computing, 08 2000.

[4] M. Conder. Trivalent (cubic) symmetric graphs on up to 2048 vertices, 2006. https://www.math.auckland.ac.nz/~conder/symmcubic2048list.txt

[5] S. Wilson and P. Potočnik. Recipes for edge-transitive tetravalent graphs. The
Art of Discrete and Applied Mathematics, 3, 08 2016.

[6] S. Wilson and P. Potočnik. Census of edge-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/TetraSS2016/TetraSS2016.zip 

[7] P. Potočnik, P. Spiga, and G. Verret. Census of arc-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/work_datoteke/Census4val-640.mgm

[8] P. Potočnik. A list of 4-valent 2-arc-transitive graphs and finite faithful amalgams of index (4, 2). European Journal of Combinatorics, 30(5):1323–1336, 2009. Part Special Issue on Metric Graph Theory.

[9] P. Potočnik, P. Spiga, and G. Verret. Census of 2-arc-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/research/Census4val2AT-2000.mgm

[10] P. Potočnik. Current list for “census of pentavalent arc-transitive graphs”. URL: https://users.fmf.uni-lj.si/potocnik/work_datoteke/AT5-Census.mgm

