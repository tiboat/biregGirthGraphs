# Computational methods for finding bi-regular graphs of given girth

This repository contains code related to finding $(\lbrace r,m\rbrace ;g)$-graphs, which are graphs of girth $g$ with vertices of degree $r$ and $m$, where $r < m$. This code is related to "J. Goedgebeur, J. Jooken and T. Van den Eede. Computational methods for finding bi-regular cages, manuscript."

There are three methods provided in the following subdirectories.

 - `biregGirthGen`: exhaustively generate all $(\lbrace r,m\rbrace ;g)$-graphs of order $n$, given $r$, $m$, $g$ and $n$.
 - `regToBireg`: try to construct $(\lbrace r,m\rbrace ;g)$-graphs starting from $r$-regular graphs.
  - `calcMaxAtMinDist`: calculate the maximum amount of vertices or edges at miminimum pairwise distance $\left\lceil g/2\right\rceil$ in an $(r,g)$-graph ($=$ $r$-regular graph of girth $g$).

For further details on these three methods, we refer to the README in the corresponding subdirectory.


## Resources

Parts of the code are taken from or based on the following resources.
 - [GenK2Hypohamiltonian](https://github.com/JarneRenders/GenK2Hypohamiltonian): the $\texttt{bitset}$ implementation, $\texttt{bitset}$ related functions and datastructures, code for converting to a graph to [graph6](http://users.cecs.anu.edu.au/~bdm/data/formats.txt) format and smaller code parts are based on this.
 - [edgeGirthRegularGraphs](https://github.com/JorikJooken/edgeGirthRegularGraphs): for an implementation of updating distances in `biregGirthGen`.
 - [cbitset](https://github.com/lemire/cbitset/tree/master): a dynamically allocated bitset implementation used in `regToBireg`.
 - [nauty](https://pallini.di.uniroma1.it/): for checking graph isomorphisms and other functionality regarding the manipulation of graphs and computing graph properties. Version 2.8.8 of nauty is included in the subdirectory `nauty`.

## Author

Tibo Van den Eede

E-mail: tibo [dot] vandeneede [at] kuleuven [dot] be
