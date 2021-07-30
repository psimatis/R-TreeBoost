# rTreeBoost
An example of using the Boost library's RTree.

CAPACITY is the max capacity of the nodes.

Creates the index from a data file.

Executes both range and knn queries.

Size is calculated using a custom visitor. The visitor is inspired by here: http://boost-geometry.203548.n3.nabble.com/How-could-I-get-nodes-MBRs-of-the-R-Tree-td4026812.html

In order to count the nodes visited by the queries, 2 files in Boost are modified. I include the files here (distance_query.hpp and spatial_query.hpp) but its possible it won't work in future Boost versions. Boost version used for this code is 1.75.0.
