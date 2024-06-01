# Algorithms for finding biregular graphs with given girth

This repository contains code related to finding $(\{r,m\};g)$-graphs, which are graphs of girth $g$ with vertices of degree $r$ and $m$.

There are two methods provided in the following subdirectories.

 - [biregCageGen](biregCageGen/): exhaustively generate all $(\{r,m\};g)$-graphs of order $n$, given $r$, $m$, $g$ and $n$.
 - [regToBireg](regToBireg/): try to construct $(\{r,m\};g)$-graphs starting from $r$-regular graphs.

For further details on these two methods, we refer to the README in the corresponding subdirectory.


## Resources

Parts of the code are taken from or based on the following resources.
 - [GenK2Hypohamiltonian](https://github.com/JarneRenders/GenK2Hypohamiltonian): the $\texttt{bitset}$ implementation, $\texttt{bitset}$ related functions and datastructures, code for converting to a graph to [graph6](http://users.cecs.anu.edu.au/~bdm/data/formats.txt) format and smaller code parts are based on this.
 - [edgeGirthRegularGraphs](https://github.com/JorikJooken/edgeGirthRegularGraphs): for an implementation of updating distances in [biregCageGen](biregCageGen/).
 - [cbitset](https://github.com/lemire/cbitset/tree/master): a dynamically allocated bitset implementation used in [regToBireg](regToBireg/).
 - [nauty](https://pallini.di.uniroma1.it/): for checking graph isomorphisms and other functionality regarding the manipulation of graphs and computing graph properties. Version 2.8.8 of nauty is included in the subdirectory `nauty`.

## Author

Tibo Van den Eede

E-mail: tibo [dot] vandeneede [at] student [dot] kuleuven [dot] be