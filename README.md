# sRanges

- [Introduction](#introduction)
- [Data](#data)
- [Output](#output)
- [Usage](#usage)
- [Tree Summarization](#tree-summarization)
- [Model Features](#model-features)
  - [Tree Encoding](#tree-encoding)
  - [Clock and Substitution Models](#clock-and-substitution-models)
- [Contributing](#contributing)
- [License](#license)

## Introduction

sRanges is a BEAST2.7 package for performing phylodynamic inference under the birth-death stratigraphic ranges model. Depending on data interpretation it can be used for both macroevolution and epidemiology.

THe dependencies on other BEAST2.7 packages can be found in `version.xml` file.

## Data

sRanges uses three types of data to inform the inference:

- **Genetic sequence alignment**
- **Morphological characters alignment**
- **Sampling data, including range of a taxa**

Sampling at present and through time is allowed, as well as having sampled ancestors. The data includes the time and uncertainty of these samples. Stratigraphic ranges are represented by two nodes marking the beginning and end of taxa (species or host), and data and time uncertainty can be associated with either of these nodes. Samples of the same taxa in between are integrated out by the model and are not needed.

## Output

sRanges produces a posterior distribution of model parameters and oriented trees:

- In the macroevolution case, they are interpreted as a budding speciation model, where each child of a node is either an ancestor or descendant species.
- In epidemiology, they are transmission trees where each child of a node is either a donor or recipient.

## Usage

You can install the latest release by adding the link https://raw.githubusercontent.com/jugne/stratigraphic-ranges/master/package.xml as a third party BEAST package repository in Beauti and installing the sRanges package that appears. 

The Beauti template is currently not available. You can find examples of BEAST2 XML in the `examples` folder.

## Tree Summarization

The package includes **SR Tree Annotator**, a tool for summarizing posterior tree samples using relationship-based credibility instead of traditional clade-based methods. This preserves the orientation information unique to SR trees.

### Launching from App Launcher (GUI)

1. Open BEAST2 App Launcher
2. Select **SR Tree Annotator** from the list
3. Fill in the input/output files and options via the GUI

### Command Line Usage

```bash
java sr.treeannotator.SRTreeAnnotator -trees input.trees -out output.tree -burnin 10
```

### Options

| Option | Description |
|--------|-------------|
| `-trees <file>` | Input tree file (Nexus or Newick format) |
| `-out <file>` | Output file for the MCC tree |
| `-burnin <n>` | Burnin percentage (default: 10) |
| `-sumProbabilities true` | Use sum instead of log product for scoring |
| `-summary <file>` | Write relationship summary to file |
| `-detailed true` | Include detailed relationship annotations |

### Example

```bash
java sr.treeannotator.SRTreeAnnotator \
    -trees posterior.trees \
    -out mcc.tree \
    -burnin 20 \
    -summary relationships.txt \
    -detailed true
```

### Output Annotations

The MCC tree is annotated with:
- `posterior` - relationship probability at each node
- `height_mean`, `height_median`, `height_95%_HPD` - height statistics

With `-detailed true`, additional annotations are added:

**For ancestry relationships (sampled ancestor nodes):**
- `relationship_type` = "ancestry"
- `ancestor_taxon` - the ancestor taxon name
- `descendant_taxa` - set of descendant taxa, e.g., `{A,B,C}`

**For orientation relationships (bifurcation nodes):**
- `relationship_type` = "orientation"
- `ancestral_taxa` - set of taxa in the ancestral lineage
- `descendant_taxa` - set of taxa in the descendant lineage

### Summary File Output

With `-summary <file>`, a text file is written containing all relationships found:

```
Ancestry Relationships:
  (A, {B,C,D}): count=850, prob=0.8500
  (E, {F,G}): count=620, prob=0.6200

Orientation Relationships:
  ({A,B} => {C,D}): count=720, prob=0.7200
  ({E} => {F,G}): count=550, prob=0.5500
```

Where `count` is the number of trees containing that relationship and `prob` is the posterior probability.

## Model Features

### Tree Encoding

We model a labeled oriented tree on stratigraphic ranges. It contains metadata on the orientation of internal nodes, where the left child represents the ancestral species (donor), and the right child represents the descendant species (recipient). For every branch, it's metadata is recorded at the child node and represents the ancestral/descendant assignment. Sampled ancestors are always attached from the left, and appropriate metadata is recorded.

Note that while trees are logged with correct node orientation, metadata should be referenced when allocating branch label.

### Clock and Substitution Models

Clock and substitution models implemented in BEAST2 are allowed.

## For developers

### Repository structure

### Contributing

Pull requests are welcome. If you have major suggestions, please discuss them before any real effort is put into development.

## Acknowledgements

The sRanges are developed by Alexandra (Sasha) Gavryushkina and Ugne Stolz. 

## License

[MIT](https://choosealicense.com/licenses/mit/)

