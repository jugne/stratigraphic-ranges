## Table of Contents

- [Tools](#tools)
- [TreeWithMetadataLogger](#treewithmetadatalogger)

# Tools 
<a name="#tools"></a>

The `Tools` class provides various utility methods for working with sRange trees.

## Class Description

The `Tools` class contains static methods that perform the following tasks:

- `getFirstLeafID(Node node)`: Retrieves the ID of the first leaf node in the tree. If the provided node is a leaf, its ID is returned. Otherwise, the method sorts the leaf nodes by height and returns the ID of the leaf node with the highest height.
- `removeLastSubstring(String split, String s)`: Removes the last substring from a string based on a specified delimiter. The method splits the string using the delimiter, removes the last substring, and then joins the remaining substrings back together using the delimiter.
- `orientateNodeChildren(int subtreeRootNr, Tree tree)`: Orientates the children of a subtree's root node based on stored metadata. If the metadata indicates that the children are not in the correct orientation, the method swaps their positions.
- `equalHeightWithPrecision(Node n1, Node n2)`: Compares the heights of two nodes with a precision threshold. If the absolute difference between the heights is within the threshold, the method returns true; otherwise, it returns false.
- `getAncestralRange(Node node)`: Retrieves the ancestral range of a node. If the node is a leaf, the method returns the node's range. Otherwise, the method returns the range of the node's first child.

# TreeWithMetadataLogger
<a name="#treewithmetadatalogger"></a>

The `TreeWithMetadataLogger` class is used to log sRange trees with associated metadata.

## Class Description

The `TreeWithMetadataLogger` class extends `BEASTObject` and implements the `Loggable` interface. It provides methods to log the sRange tree with metadata.

The class contains the following inputs:

- `srTreeInput`: The input for the sRange tree to be logged (required).
- `clockModelInput`: The input for the branch rate model to be logged with the branches of the tree.
- `substitutionsInput`: The input to specify whether to report branch lengths as substitutions.
- `decimalPlacesInput`: The input to specify the number of decimal places to use when writing branch lengths and rates.
- `parameterInput`: The input for metadata to be logged with the tree nodes.
- `logOrientationInput`: The input to specify whether to report if a node is a donor or recipient.

The `TreeWithMetadataLogger` class provides the following functionality:

- Logging the sRange tree with associated metadata.
- Formatting branch lengths and rates based on the specified decimal places.
- Building the Newick string representation of the tree with metadata.
- Appending doubles to the output buffer with appropriate formatting.
