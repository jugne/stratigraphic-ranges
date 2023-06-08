This folder contains operators for trees with stratigraphic ranges. Below the code structure and functionality for each operator is outlined.

# Table of Contents
- [SRTreeOperator](#srtreeoperator)
- [SRWilsonBalding](#srwilsonbalding)
- [SRLeafToSampledAncestorJump](#srleaftosampledancestorjump)
- [LeftRightChildSwap](#leftrightchildswap)

# SRTreeOperator

The `SRTreeOperator` class is an abstract class that extends the `Operator` class and serves as a base class for operators in the sRange tree.

## Class Description

The `SRTreeOperator` class extends the `Operator` class and provides common functionality for operators in the sRange tree. It contains an input for the sRange tree on which the operation is performed and an input for marking ancestors of changed nodes as changed.

### Methods

- `getOtherChild(parent, child)`: Returns the other child node of a given parent node.
- `replace(parent, child, replacement)`: Replaces a child node with another node in the parent node.

# SRWilsonBalding

The `SRWilsonBalding` class implements the Wilson-Balding proposal for the sRange tree.

## Class Description

The `SRWilsonBalding` class extends the `SRTreeOperator` class and implements the Wilson-Balding proposal for the sRange tree.

### Methods

- `proposal()`: Performs the Wilson-Balding proposal and returns the logarithm of the Hastings ratio.
# SRLeafToSampledAncestorJump

The `SRLeafToSampledAncestorJump` class implements a narrow move between trees of different dimensions in the sRange tree.

## Class Description

The `SRLeafToSampledAncestorJump` class extends the `SRTreeOperator` class and implements a narrow move between trees of different dimensions in the sRange tree. It selects a random leaf node with a younger sibling or a sampled internal node and performs the move by either replacing the leaf with a sampled internal node or inserting a new parent node at a height uniformly chosen between the sampled node's height and its old parent's height.

### Methods

- `proposal()`: Performs the narrow move and returns the logarithm of the Hastings ratio.

# LeftRightChildSwap

The `LeftRightChildSwap` class is an operator that swaps the left and right child nodes of a randomly chosen internal node in the sRange tree.

## Class Description

The `LeftRightChildSwap` class extends the `SRTreeOperator` class and implements the `proposal` method to perform the child swap operation. It chooses a random internal node in the tree (excluding the root and leaf nodes) that does not belong to the same sRange as its left child, and swaps its left and right child nodes.

### Methods

- `proposal()`: Performs the left-right child swap operation by swapping the left and right child nodes of a randomly chosen internal node.
