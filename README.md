# sRanges

sRanges is a BEAST2 package for performing phylodynamic inference under both birth-death stratigraphic ranges model for macroevolution and epidemiology. It produces a posterior distribution of model parameters and oriented trees:
    - in macroevolution case they are interpreted as budding speciation model, where each child of a node is either ancestor or descendant species. 
    - in epidemiology, they are transmission trees where each child of a node is either donor or recipient. 
   

## Data

sRanges use three types of data to inform the inference:

#### Genetic sequence alignment
#### Morphological characters alignment 
#### Sampling data
Sampling at present and through time is allowed, as well as having sampled ancestors. The data is the time and uncertainty of these samples. 

##### Stratigraphic ranges
The range of taxa (species or host) is represented by two nodes marking its beginning and end. Data and time uncertainty can be associated with either of these nodes. Same taxa samples in between are integrated out by the model and are not needed.    



## Output


## Usage

The Beauti template is currently not available. You can find examples in the `examples` folder. 


## Model features

### Tree encoding

### Clock and substitution models models

Clock and substitution models which are implemented in BEAST2 are allowed.


## Contributing

Pull requests are welcome. If you have major suggestions, please discuss it before any real effort is put into development.

## License

[MIT](https://choosealicense.com/licenses/mit/)
