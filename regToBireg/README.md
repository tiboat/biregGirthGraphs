# regToBireg

This directory contains code to try to construct $(\lbrace r,m\rbrace ;g)$-graphs based on $r$-regular graphs. Given a list of $r$-regular graphs, $m$, $g$ and the range $[minN, maxN]$ in which the order of the $(\lbrace r,m\rbrace ;g)$-graph should be, the code loops over the list. For each $r$-regular graph $G$ it first checks three necessary conditions on the order, girth and diameter of $G$, to make it feasible to obtain an $(\lbrace r,m\rbrace ;g)$-graph with order $n \in [minN, maxN]$. Only if these three conditions are met, then the construction is tried. There are three constructions, that each have a different condition on the value on $m$.

 - **Construction 1**: $m=r+1$. The code is available in `regToBiregAddEdge.c`.
 - **Construction 2**: $m$ is even. The code is available in `regToBireg.c`.
 - **Construction 3**: $m$ is odd. The code is available in `regToBireg.c`.
## Data

### Lists of regular graphs

This directory contains compressed versions of lists or regular graphs. More precisely, 
 - Files starting with `cubic_arc_trans` contain $3$-regular arc-transitive graphs, where the two following numbers indicate the range of orders of the graphs in the file [1]. For example `cubic_arc_trans_4-4886.zip` contains $3$-regular arc-transitive graphs of order $n \in [4,4886]$.
 - `cubic_semi_symm-10000.zip` contains $3$-regular semi-symmetric graphs up to 10000 vertices [1]. 
 -  `quadric.zip` contains lists of $4$-regular graphs which are edge-transitive, arc-transitive or 2-arc-transitive [6,7,8,10,11].
 -  `penta_arc_trans.zip` contains a list of $5$-regular arc-transitive graphs [5].
   
Unzip these files in order to use the lists inside of them. Also, other lists where used in [3]. These are the following.

 - The repository [cubic-vertextransitive-graphs](https://github.com/kguo-sagecode/cubic-vertextransitive-graphs) contains $3$-regular vertex-transitive graphs in sparse6 format up to 1280 vertices. This data originates from [4,9].
 - On [this link](https://users.fmf.uni-lj.si/potocnik/CubicCay/CubicCayUpTo4094.zip) (accessible via [graphsym.net](https://graphsym.net) ) a .zip file (7.2 GB) of $3$-regular Cayley graphs up to 4094 vertices can be downloaded in sparse6 format [4].

### Obtained record graphs

The directory `found_graphs` contains $(\lbrace r,m\rbrace ;g)$-graphs encoded in sparse6 format of smaller order than we found in the literature and that did not lead to a better improvement than Theorem 4.2 of [3]. Each graph is stored in a different file formatted as `r<r>_m<m>_g<g>_n<n>.s6`, e.g. `r3_m4_g14_n544.s6`.


## Code

### Compilation


First configure and compile nauty. Go to the nauty directory (`biregGirthGraphs/nauty/`) and execute the following two commands.

```
./configure
make
```

Then compile `regToBiregAddEdge.c` and `regToBireg.c` by using the `make` command, which creates the executables `regToBiregAddEdgeExec` and `regToBiregExec`.
```
make
```

### Execution

Then one should feed the executable with a list of $r$-regular graphs via stdin and give the executable the proper arguments.
```
cat <list_file> | ./regToBiregAddEdgeExec <r> <g> <minN> <maxN>
````
```
cat <list_file> | ./regToBiregExec <r> <m> <g> <minN> <maxN>
````
where `<list_file>` is a list of $r$-regular graphs in sparse6 or graph6 format and `<r>`, `<m>`, `<g>`, `<minN>` and `<maxN>` have the same meaning as described above. Note that we don't need to specify `<m>` for `regToBiregAddEdgeExec` since there $m=r+1$.

#### Example
```
cat graph_lists/quadric-arc-trans.g6 | ./regToBiregExec 4 7 9 281 742
````
This goes through the list of $3$-regular vertex-transitive graphs of order in $[281,742]$ in order to try to construct a $(\lbrace 4,7\rbrace ;9)$-graph of order $n \in [281, 742]$.

### Output
The code outputs the order, girth and diameter of the regular graphs that satisfy the three conditions. If an $(\lbrace r,m\rbrace ;g)$-graph of order $n \in [minN, maxN]$ is found by the construction, this graph is printed in sparse6 format and the code does not continue trying to construct more graphs.

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
:~?Ee_K?gD?Q [...] Ph_uc
filtered 2586 graphs out of 2590
```
where the complete sparse6 string is trimmed to `:~?Ee_K?gD?Q [...] Ph_uc`. 



## References

[1] M. Conder. Trivalent symmetric graphs on up to 768 vertices. Journal of
Combinatorial Mathematics and Combinatorial Computing, 08 2000.

[2] M. Conder. Trivalent (cubic) symmetric graphs on up to 2048 vertices, 2006. https://www.math.auckland.ac.nz/~conder/symmcubic2048list.txt

[3] J. Goedgebeur, J. Jooken and T. Van den Eede. Computational methods for finding bi-regular cages, manuscript.

[4]     P. Potočnik, P. Spiga, G. Verret, Cubic vertex-transitive graphs on up to 1280 vertices, Journal of Symbolic Computation 50 (2013), 465-477. 

[5] P. Potočnik. Current list for “Census of pentavalent arc-transitive graphs”. URL: https://users.fmf.uni-lj.si/potocnik/work_datoteke/AT5-Census.mgm

[6] P. Potočnik. A list of 4-valent 2-arc-transitive graphs and finite faithful amalgams of index (4, 2). European Journal of Combinatorics, 30(5):1323–1336, 2009. Part Special Issue on Metric Graph Theory.

[7] P. Potočnik, P. Spiga, and G. Verret. Census of 2-arc-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/research/Census4val2AT-2000.mgm

[8] P. Potočnik, P. Spiga, and G. Verret. Census of arc-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/work_datoteke/Census4val-640.mgm

[9] P. Potočnik, P. Spiga, and G. Verret. Bounding the order of the vertex-stabiliser in 3-valent vertex-transitive and 4-valent arc-transitive graphs. Journal of Combinatorial Theory, Series B, 111:148–180, 2015.

[10] S. Wilson and P. Potočnik. Census of edge-transitive tetravalent graphs. https://users.fmf.uni-lj.si/potocnik/TetraSS2016/TetraSS2016.zip 

[11] S. Wilson and P. Potočnik. Recipes for edge-transitive tetravalent graphs. The Art of Discrete and Applied Mathematics, 3:41–63, 2016.













