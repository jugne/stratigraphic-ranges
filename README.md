# sRanges

- [Introduction](#introduction)
- [Data](#data)
- [Output](#output)
- [Usage](#usage)
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

