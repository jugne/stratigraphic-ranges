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

sRanges is a BEAST2 package for performing phylodynamic inference under both the birth-death stratigraphic ranges model for macroevolution and epidemiology.

## Data

sRanges uses three types of data to inform the inference:

- **Genetic sequence alignment**
- **Morphological characters alignment**
- **Sampling data**

Sampling at present and through time is allowed, as well as having sampled ancestors. The data includes the time and uncertainty of these samples. Stratigraphic ranges are represented by two nodes marking the beginning and end of taxa (species or host), and data and time uncertainty can be associated with either of these nodes. Samples of the same taxa in between are integrated out by the model and are not needed.

## Output

sRanges produces a posterior distribution of model parameters and oriented trees:

- In the macroevolution case, they are interpreted as a budding speciation model, where each child of a node is either an ancestor or descendant species.
- In epidemiology, they are transmission trees where each child of a node is either a donor or recipient.

## Usage

The Beauti template is currently not available. You can find examples of BEAST2 XML in the `examples` folder.

## Model Features

### Tree Encoding

### Clock and Substitution Models

Clock and substitution models implemented in BEAST2 are allowed.

## Contributing

Pull requests are welcome. If you have major suggestions, please discuss them before any real effort is put into development.

## License

[MIT](https://choosealicense.com/licenses/mit/)

