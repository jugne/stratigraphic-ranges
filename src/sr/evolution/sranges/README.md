# StratigraphicRange

The `StratigraphicRange` class represents a range of occurrences for a taxon in a phylogenetic tree.

## Table of Contents
- [Package](#package)
- [Class Description](#class-description)
- [Input Fields](#input-fields)
- [Instance Variables](#instance-variables)
- [Methods](#methods)

## Package

`sr.evolution.sranges`

## Class Description

The `StratigraphicRange` class extends the `BEASTObject` class and is used to define the first and last occurrences of a taxon and maintain a list of nodes representing the range.

### Input Fields

- `taxonFirstOccurrenceInput`: Represents the first occurrence of the taxon as a BEAST `Taxon` object.
- `taxonLastOccurrenceInput`: Represents the last occurrence of the taxon as a BEAST `Taxon` object.

### Instance Variables

- `firstOccurrenceID`: Stores the ID of the first occurrence taxon.
- `lastOccurrenceID`: Stores the ID of the last occurrence taxon.
- `isSingleFossilRange`: Indicates whether the range represents a single fossil occurrence.
- `nodes`: A list of integers representing the nodes that belong to the range.
- `directAncestorNodes`: A list of integers representing the direct ancestor nodes in the range.

### Methods

- `initAndValidate()`: Initializes and validates the object, checking if both the first and last occurrence taxa are specified.
- `containsNodeNr()`: Checks if a given node number is present in the range.
- `addNodeNrAfter()`: Adds a node number after a specified node in the range.
- `removeNodeNr()`: Removes a node number from the range.
- `removeAllNodeNrs()`: Clears all node numbers from the range.
- `getNodeNrs()`: Retrieves the list of node numbers in the range.
- `setFirstOccurrenceNodeNr()`: Sets the node corresponding to the first occurrence of the range.
- `setLastOccurrenceNodeNr()`: Sets the node corresponding to the last occurrence of the range.
- `makeSingleFossilRange()`: Converts the range to represent a single fossil occurrence.
- `isSingleFossilRange()`: Checks if the range represents a single fossil occurrence.
- `getInternalNodeNrs()`: Retrieves the node numbers of internal nodes in the range.
- `setFirstOccurrenceID()`: Sets the ID of the first occurrence taxon.
- `setLastOccurrenceID()`: Sets the ID of the last occurrence taxon.
- `getFirstOccurrenceID()`: Retrieves the ID of the first occurrence taxon.
- `getLastOccurrenceID()`: Retrieves the ID of the last occurrence taxon.
- `addNodeNr()`: Adds a node number to the range.
