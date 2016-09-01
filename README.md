# PyPedWorks (or PedWorks.py)
A Python program which applies network theory to model genealogical structured data (pedigrees).

**Authors:**

* Bruno C. Perez <brunocpvet@gmail.com>
* Ricardo V. Ventura
* Julio C. C. Balieiro

**Latest version:** 0.1.9 - (Sep 1, 2016)

**Requires:**

* NetworkX
* Numpy
* Pandas
* Collections
* Louvain [download](http://perso.crans.org/aynaud/communities/)


### Instructions:
  * Open the terminal.

  * Go to pedWorks.py file directory and run:


``python2 pedWorks.py test.ini`` - for systems where Python 3.X is default.

or

``python pedWorks.py test.ini`` - for systems where Python 2.7 is default.

## Functions:

### Pedigree reordering/renumbering by Kahn's algorithm

### Network analysis

Performs a network analysis after applying network theory into a pedigree.

* Calculates basic measures: Number of nodes/edges, Density, Average (out)degree
* Calculates centralities for all nodes:
       * (out)Degree centrality
       * Closeness Centrality
       * Betweeness Centrality
       * Eigenvector Centrality
       * If ``ctplot = TRUE``, creates a plot for each centrality.
* Creates the (out) degree rank plot for the pedigree.
    

### Pedigree (as network) drawing

Uses Fruchterman & Reingold directed force drawing algorithm for a pedigree modeled as a acyclic directed network.

* ``draw_simple`` - Draws a simple one-colored network.
* ``draw_group`` - Draws distinct groups of animals in different colors.
* ``draw_multigroup`` - Draws only pre-ditermined nodes (individuals) from the pedigree maintainning their relationship as calculated for the complete pedigree.
* ``draw_cluster`` - Uses Louvain method for community detection to find the underlying structure of a pedigree.

### Breed Composition calculation

Makes use of a network framework to calculate the breed composition of all animals in a multi-breed population pedigree.


## Basic .ini file options:
A brief description of the available options.

  * ``Reorder`` function arguments:
     * ``outfile`` = specifies the name of the reordered pedigree file. ex: ``pedOrd.txt``
     * ``outheader`` = ``TRUE/FALSE`` if ``TRUE``, prints a header on the ordered pedigree file.
     * ``format`` = specifies the format for the ordered pedigree output.
       * ``fwf`` - fixed width format
       * ``csv`` - comma separated format
       * ``txt`` - space separated format


  * ``Analysis`` function arguments:
     * ``ctplot`` = ``TRUE/FALSE`` if ``TRUE``, creates the centrality plots for the pedigree.
       * (Out)Degree, Closeness, Betweeness and Eigenvector centralities.
 
 
  * ``Draw`` function arguments:
     * ``nscale`` = specifies scale in which the network will be plotted (the x and y axis range)
     * ``nalpha`` = specifies the transparency of the nodes. ex:. ``0.1`` - ``1.0``
     * ``nsize`` = specifies the size of the nodes. ex:. ``10`` - ``100``.
     * ``ncolor`` = specifies the RGB code for the color of the nodes. [RBG colochart](http://www.rapidtables.com/web/color/RGB_Color.htm)
     * ``ealpha`` = specifies the transparency of the edges. ex:. ``0.1`` - ``1.0``
     * ``ewidth`` = specifies the width of the edges. ex:. ``0.1`` - ``1.0``
     * ``ecolor`` = specifies the RGB code for the color of the nodes. [RBG colochart](http://www.rapidtables.com/web/color/RGB_Color.htm)
     * ``initpos`` = ``TRUE/FALSE`` - If ``TRUE``, node positioning start as defined in a ``.txt`` file.
     * ``posfile`` = specifies the file (``.txt``) containning node positioning. (format: ``node  x-axis  y-axis``)
     * ``savepos`` = ``TRUE/FALSE`` - If True, saves final node positioning in ``nodepos.txt`` file.

* ``draw_ group`` and ``draw_multigroup`` function specific argument:
     * ``group_list`` = specifies he name(s) of the file(s) containning a list of individuals to be highlighted. It may received multiple groups. ex:. ``[group1.txt, group2.txt, group3.txt]``
     * ``color_list`` = specifies the RGB code for the color for each highlighted group. May receive multiple colors that must be specified in the same number and order of the groups in ``group_list``. ex:. ``[#008000, #4682B4, #DC143C]``

* ``draw_cluster`` function specific argument:
     * ``cSize`` = specifies the community minimum size treshold. Communities containning less than cSize animals will not be colored in the network drawing. ex: ``5`` - ``20``

  * ``Breed Composition`` function arguments:
     * ``infile`` = specifies the input pedigree file. ex: ``pedigree.txt``
     * ``nbreed`` = specifies the number of different breeds in the population. ex: ``2`` - ``4``.


## Pedigree Network drawing examples
 The following examples were obtained by real genealogical data analysis using PedWorks. 
  * Topological sorting by Fruchterman & Reingold force-directed method.
  * Community detection algorithm by Louvain method.
  * Nodes (individuals) with same colors belong to the same community.
 
### Example Pedigree [1]
 * 11,767 individuals (cattle)
 * Mean inbreeding coefficient of 2.30%
 * Number of distinct relevant communitities = 31
[![example1.png](https://s14.postimg.org/v2w1gamqp/example1.png)](https://postimg.org/image/ibhv9scyl/)

### Example Pedigree [2]
 * 8,425 individuals (cattle)
 * Mean inbreeding coefficient of 1.60%
 * Number of distinct relevant communitities = 20
[![example2.png](https://s9.postimg.org/5va1e1jb3/example2.png)](https://postimg.org/image/r4xnovzln/)
