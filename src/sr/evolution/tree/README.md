This folder contains classes for sRanges tree and nodes. Below the code structure and functionality for each class is outlined.



# Table of Contents
- [SRTree](#srtree)
- [SRNode](#srnode)
- [RandomSRangeTree](#randomsrangetree)
- [SampledAncestorLogger](#sampledancestorlogger)
- [SpeciationLogger](#speciationlogger)

# SRTree

The `SRTree` class represents a labeled oriented tree on stratigraphic ranges under budding speciation. It extends the `Tree` class and implements the `TreeInterface`.

## Class Description

The `SRTree` class represents a labeled oriented tree on stratigraphic ranges. It contains metadata on the orientation of internal nodes, where the left child represents the ancestral species, and the right child represents the descendant species. Sampled ancestors are always attached from the left, and appropriate metadata is recorded.

### Input Fields

- `stratigraphicRangeInput`: Represents all stratigraphic ranges as a list of `StratigraphicRange` objects.
- `treeInput`: Represents the tree to start with as a `Tree` object.

### Instance Variables

- `sRanges`: An `ArrayList` of `StratigraphicRange` objects representing the stratigraphic ranges.
- `storedSRanges`: An `ArrayList` of stored `StratigraphicRange` objects.

### Methods

- `initAndValidate()`: Initializes and validates the object, assigns the tree if provided, and initializes the stratigraphic ranges.
- `initSRanges()`: Initializes the stratigraphic ranges based on the input or inferred from the provided tree structure.
- `initStoredRanges()`: Initializes the stored stratigraphic ranges based on the current stratigraphic ranges.
- `assignFrom()`: Copies all values from an existing tree.
- `assignFromFragile()`: Copies the tree structure only.
- `assignFrom()`, `assignFromFragile()`, `assignFrom()` (helper): Helper methods for tree assignment.
- `store()`: Stores the current state of the tree and stratigraphic ranges.
- `restore()`: Restores the tree and stratigraphic ranges from the stored state.
- `fromXML()`: Populates the tree from XML file data.
- `orientateTree()`: Orients the tree according to stored metadata.
- `addOrientationMetadata()`, `addOrientationMetadataNode()`: Adds orientation metadata to each node.
- `getSRanges()`: Retrieves the stratigraphic ranges of the tree.
- `getNonSingleSRangeCount()`: Retrieves the count of non-single fossil stratigraphic ranges.
- `getSRangesInternalNodeNrs()`: Retrieves the internal node numbers in the stratigraphic ranges.
- `sRangesContainsID()`: Checks if the stratigraphic ranges contain a specific taxon ID and output the range containing it.
- `getRangeOfNode()`: Retrieves the stratigraphic range to which a node belongs.
- `getSharedRange()`: Retrieves the shared stratigraphic range between two nodes.
- `belongToSameSRange()`: Checks if two nodes belong to the same stratigraphic range.
- `log()`: Logs the state of the tree.

# SRNode

The `SRNode` class represents a node in a labeled oriented tree on stratigraphic ranges. It extends the `Node` class.

## Class Description

The `SRNode` class represents a node in a labeled oriented tree on stratigraphic ranges.

### Methods

- `copy()`: Returns a deep copy of the node.
- `sort()`: Throws a runtime exception indicating that ordered trees should not be sorted.
- `toShortNewickForLog(boolean printInternalNodeNumbers)`: Returns the tree in Newick format with length and metadata information. Nodes are numbered instead of using node labels. Internal nodes with non-null IDs are also printed. If `printInternalNodeNumbers` is set to true, all internal nodes are labeled, which is useful for storing the state to a file.

# RandomSRangeTree

The `RandomSRangeTree` class extends the `SRTree` class and implements the `StateNodeInitialiser` interface. It represents a random tree with stratigraphic ranges.

## Class Description

The `RandomSRangeTree` class represents a random tree with stratigraphic ranges.

### Inputs

- `taxaInput`: Input of type `Alignment` representing the set of taxa to initialize the tree specified by the alignment.
- `populationFunctionInput`: Input of type `PopulationFunction` representing the population function for generating coalescent (required).

### Fields

- `nrOfTaxa`: Total number of taxa.
- `children`: An array of lists to store the children of each node.
- `taxa`: A set to store the taxa.
- `nextNodeNr`: Number of the next internal node, used when creating new internal nodes.

### Methods

- `initAndValidate()`: Initializes and validates the tree by setting up taxa, initializing state nodes, and calling the `initAndValidate()` method of the superclass.
- `swap(List list, int i, int j)`: Swaps elements in a list at indices `i` and `j`.
- `initStateNodes()`: Initializes state nodes by setting up taxa, simulating the tree, and assigning node numbers.
- `setNodesNrs(Node node, int internalNodeCount, int[] n, Map<String, Integer> initial)`: Sets the numbers of nodes recursively.
- `getInitialisedStateNodes(List<StateNode> stateNodes)`: Gets the initialized state nodes.
- `simulateTree(Set<String> taxa, PopulationFunction demoFunction)`: Simulates a coalescent tree given a taxon list and a demographic function.
- `processCandidateTraits(Set<Node> candidates, List<TraitSet> traitSets)`: Applies traits to a set of nodes.
- `simulateCoalescent(Set<Node> candidates, PopulationFunction demoFunction)`: Simulates coalescent events given a set of candidates and a demographic function.
- `simulateCoalescent(List<Node> nodes, PopulationFunction demographic, double currentHeight)`: Simulates coalescent events given a list of nodes, a demographic function, and a current height.
- `getMinimumInactiveHeight()`: Returns the height of the youngest inactive node.
- `setCurrentHeight(double height)`: Sets the current height.
- `getActiveNodeCount()`: Returns the number of active nodes.
- `coalesceTwoActiveNodes(double height)`: Coalesces two nodes in the active list by removing them and replacing them with a new node.
- `nodeList`: An ArrayList to store the nodes.
- `activeNodeCount`: The number of active nodes.

# SampledAncestorLogger

The `SampledAncestorLogger` class extends the `CalculationNode` class and implements the `Loggable` and `Function` interfaces. It represents a logger for reporting sampled ancestor (SA) count in a tree.

## Class Description

The `SampledAncestorLogger` class is responsible for logging the SA count in a tree.

### Inputs

- `treeInput`: Input of type `SRTree` representing the tree to report the SA count for (required).

### Methods

- `initAndValidate()`: Initializes and validates the logger.
- `init(PrintStream out)`: Initializes the logger for logging.
- `log(long nSample, PrintStream out)`: Logs the SA count.
- `close(PrintStream out)`: Closes the logger.
- `getDimension()`: Returns the dimension of the logged data (1 in this case).
- `getArrayValue()`: Returns the SA count as an array value.
- `getArrayValue(int iDim)`: Returns the SA count as an array value for a specific dimension.

# SpeciationLogger

The `SpeciationLogger` class extends the `CalculationNode` class and implements the `Loggable` and `Function` interfaces. It represents a logger for reporting speciation events in a tree.

## Class Description

The `SpeciationLogger` class is responsible for logging speciation events in a tree.

### Inputs

- `treeInput`: Input of type `SRTree` representing the tree to report the speciation events for (required).
- `sepStringInput`: Input of type `String` representing the separator string for ranges (default: "_").
- `directStringInput`: Input of type `String` representing the string to indicate speciation direction (default: ">").
- `onlyFirstInput`: Input of type `Boolean` indicating if only the first descendant should be logged (default: true).

### Methods

- `initAndValidate()`: Initializes and validates the logger.
- `init(PrintStream out)`: Initializes the logger for logging.
- `log(long nSample, PrintStream out)`: Logs the speciation events.
- `close(PrintStream out)`: Closes the logger.
- `getDimension()`: Returns the dimension of the logged data (1 in this case).
- `getArrayValue()`: Returns the speciation count as an array value.
- `getArrayValue(int iDim)`: Returns the speciation count as an array value for a specific dimension.
