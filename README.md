# PySanger Installation and User Manual

This repository has been adapted from [ponnhide/PySanger](https://github.com/ponnhide/PySanger) to allow alignment of Sanger Sequencing files to vector maps with file endings .gb .dna or .xdna

## Installation
1.  Install the following Python packages by  
	
	```sh
	conda create -n PySanger
	conda activate PySanger
	conda config --add channels defaults
	conda config --add channels bioconda
	conda config --add channels conda-forge
	
	conda install matplotlib seaborn numpy biopython pandas logomaker snapgene-reader
	```

	```sh
	git clone https://github.com/philipptrepte/PySanger.git
	```

2.  Set PYTHONPATH to the directory where you cloned the repository.

## API
- **abi_to_dict**_`(filename=str)`_  
	Generate confidence value array and signal intensity arrays of each channel (A, T, G or C) at peak positions. 
	
	**Parameter**
	
	- **filename**: *`str`*  (default: None)  
	A file path of sanger sequencing result.   
	
	**Return**
	_dict_
	``` 
	{"conf": quality scores at peak positions.
	 "channel":{"A": signal intensities at peak positions in the channel for 'A',
	            "T": signal intensities at peak positions in the channel for 'T',
	            "G": signal intensities at peak positions in the channel for 'G',
	            "C": signal intensities at peak positions in the channel for 'C'}
	}
	```

- **generate_consensusseq**_`(abidata=dict)`_
	Generate the most consensus seq from a senger sequencing result.  

	**Parameter**
	
	- **abidata**: *`dict`*  (default: None)  
	A dict object returned by 'abi_to_dict'.   
	
	**Return**
	_tuple_
	`(str:Forward strand sequence (5'->3'), str:Reverse strand sequence (5'->3'))` 

- **generate_pwm**_`(abidata=dict)`_
	Generate position weight matrix based on signal intensities of each channel.   
	
	**Parameter**
	
	- **abidata**: *`dict`*  (default: None)  
	A dict object returned by 'abi_to_dict'.   

	**Return**
	_pandas.DataFrame_

- **alignment**_`(abidata=dict), template=str, strand=int`_
	Run alignment of sanger sequencing result against a user-provided sequence or against a vector map.

	**Parameter**  	
	
	- **abidata**: *`dict`*  (default: None)  
	A dict object returned by 'abi_to_dict'.   
	- **template**: *`str`*  (default: None)  
	`template` can be nucleotide sequence or a file path to a vector map. Allowed file extensions are `.gb`, `.dna` and `.xdna`
	- **strand**: *`str`* (1 or -1, default: 1) 
	A sequencing strand used for the alignment and visualization. `1` indicates the plus strand. `-1` indicates the minus strand. 

	**Return**
	_tuple_

- **adjust_quality_scores**_`(tasubject=str, quality_scores=_abidata)`_
	Inserts `-1` in quality scores at insertion positions observed after alignment of `abidata` with `template`

	**Parameters**

	- **tasubject**: *`str`* (default: tasubject)
	A trimmed an sequence aligned nucleotide string.

	- **quality_scores**: *`list`* (default: _abidata["conf"])
	A list of quality scores at peak positions from the sanger sequencing result.

	**Return**
	_list_

- **visualize**_`(alignment=tuple, strand=int, fig=matplolib.pyplot.figure)`_ 
	Visualize a sanger sequencing result. 
	
	**Parameter**  	
	
	- **region**: *`str`* ("all" or "aligned", default: "aligned")
	A region used for the visualization. If `all`, it will visualize the entire region of the `template` vector map and the `abidata` sequencing result. If `aligned`, it will visualize the minimum aligned region between the `template` and `abidata` sequencing result.
	
	**Return**
	_matplotlib.figure.Figure_


## Example usage 
Visualise peack intensities from a Sanger sequencing result. You can use the example by `ponnhide/PySanger` where you use as `abidate` the `.ab` file `BE MAFB5.ab1` and as `query` the sequence `AGCCGGCTGGCTGCAGGCGT`.

Alternatively, can you as `abidata` the file `seq_results/QPSQ0664-CMV-for.ab1` and as `template` the vector map file `templates/QPPL0052_pcDNA3.1_mCitrine-C1-GW.dna`

### region = "aligned"
```python
from pysanger import * 
align        = alignment(abidata="BE MAFB5.ab1", template="AGCCGGCTGGCTGCAGGCGT")
fig	         = visualize(align, region = "aligned", fontsize = 5)

fig.savefig("test_aligned.pdf", bbox_inches="tight") 
```

![test_aligned.pdf](test_aligned.pdf)


### region = "all"
``` python
from pysanger import * 
align        = alignment(abidata="seq_results/QPSQ0664-CMV-for.ab1", template="templates/QPPL0052_pcDNA3.1_mCitrine-C1-GW.dna")
fig          = visualize(align, region="all")

fig.savefig("test_all.pdf", bbox_inches="tight") 
```

![test_all.pdf](test_all.pdf)

