# rTreeBoost
An example of using the Boost library's RTree.

CAPACITY is the max capacity of the nodes.

Creates the index from a data file.

Executes both range and knn queries.

Size is calculated using a custom visitor. The visitor is inspired by here: http://boost-geometry.203548.n3.nabble.com/How-could-I-get-nodes-MBRs-of-the-R-Tree-td4026812.html
