# sRanges

sRanges is a BEAST2 package for performing phylodynamic inference under both birth-death stratigraphic ranges model for macroevolution and epidemiology. It produces oriented trees. In macroevolution case they are interpreted as budding speciation model, where each child of a node is either ancestor or descendant species. In epidemiology, there are transmission trees where each child of a node is either donor or recipient. 
    

## Data

sRanges use three types of data to inform the inference:

#### Genetic sequence alignment
#### Morphological characters alignment 
#### Sampling data
Sampling at present and through time is allowed, as well as having sampled ancestors. The data is the time and uncertainty of these samples. 

#### Stratigraphic ranges
The range of taxa (species or host) is represented by two nodes marking its beginning and end. Data and time uncertainty can be associated with both these nodes. Same taxa samples in between are integrated out by the model.    



## Output


## Usage



## Model features

### Tree encoding

### Clock models

### Substitution models




## Contributing

Pull requests are welcome. If you have major suggestions, please discuss it before any real effort is put into development.

## License

[MIT](https://choosealicense.com/licenses/mit/)