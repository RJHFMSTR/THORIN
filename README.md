# THORIN v1.2 : ```T```arget ```H```aplotype ```OR```igin ```IN```ference

## What's new in v1.2?:
THORIN v1.2 uses the same algorithm as in THORIN v1.1, but offers the opportunity to
* write the ouptut in a vcf.gz or bcf format
* directly output the IBD segments and which haplotype is maternal or paternal in these segments
* directly output a vcf with target individuals phased according to the IBD segments

To get this information, use the following parameters:
* ```--ibd``` specifies the ouptut file for the IBD segments and their states. The file contains one row per IBD segment, and its columns are structured as follow: \
**CHR | start | end | Prob | length_cM | target | UPD** \
**CHR** -> is the chromosome on which the IBD segments was found \
**start** -> is starting position of the ibd segment \
**end** -> is ending position of the ibd segment \
**Prob** -> defines which haplotype is paternal/maternal on IBD segments. It can either be A (left haplotype is paternal and right maternal), B (left haplotype is maternal and right is paternal), C (undefined) or D (Disomy). \
**length_cM** -> size of the IBD segment in cM. \
**target** -> id of the target sample.

* ```--phasing``` specifies the output file for the target phased based on their IBD sharing \

For a quick test of the new functionnalities you can run:
```
bin/thorin -I test/related.chr20.vcf.gz -H test/unrelated.chr20.vcf.gz  -G test/Trios.txt -M maps/chr20.b37.gmap.gz -R 20 -O original_ouptut.bcf --ibd ibd_output.tsv --phasing phased_target.bcf
```

## Overview

THORIN is a tool that we use to infer the PofO for alleles shared IBD with surrogate parents. It is a Hidden Markov Model (HMM) specifically design to identify IBD sharing between the target haplotypes and a reference panel mixing haplotypes from 2 different sources: from the surrogate parents of the target (labelled as mother or father) and from unrelated samples. We aimed for such a probabilistic model for its robustness to phasing and genotyping errors compared to approaches based on exact matching such as the positional Burrows–Wheeler transform (PBWT). The model then uses a forward-backward procedure to compute, for each allele of a target haplotype, the probability of copying the allele from (i) the surrogate mother haplotypes, (ii) the surrogate father haplotypes or (iii) unrelated haplotypes. We also use a set of unrelated haplotypes as decoys so that the model is not forced to systematically copy from surrogate parents. When the model copies the target haplotype from a specific surrogate parent at a given locus with high probability, we can therefore infer the PofO at this locus from the parental group the surrogate parent belongs to. When the model copies from unrelated haplotypes, no inference can be made at the locus


If you use the THORIN in your research work, please cite the following paper: [Hofmeister, R.J. et al. Parent-of-Origin inference for biobanks. Nature Communication 2022. https://doi.org/10.1038/s41467-022-34383-6](https://www.nature.com/articles/s41467-022-34383-6)


</br>
  &nbsp
  
#
## Documentation



To run an example, use the following command line, which should take less than 10 seconds:

```
./bin/thorin_static -I test/related.chr20.vcf.gz -H test/unrelated.chr20.vcf.gz -G test/Trios.txt -M maps/chr20.b37.gmap.gz -R 20 -O output_test.txt
```


</br>
  &nbsp
  
  
Mandatory Options:


```-I``` specifies the input .vcf file with all samples listed in the group file (-G).

```-H``` specifies the input .vcf file with a set of fully unrelated samples.
	
```-G``` specifies the group file.

```-M``` specifies the genetic maps.

```-O``` specifies the output file.


To see the full list of options, run the following command:
```./bin/thorin_static --help```


</br>
  &nbsp

The group file (```-G```) is formatted as follows:
* the first column specifies the target individual id
* the second column specifies the first group of relatives. It starts with ```G1=```, referring to ```Group 1```, followed by a ';'-separated list of relatives.
* the third column, if any, specifies the second group of relatives with ```G2=```.

</br>
  &nbsp

The output file stores the copying probability of each group of relatives (```G1``` and ```G2``` if any) as well as the copying probability of unrelated individuals (```H```). It is formated as follows:
* columns 1-4 specifies the chromosome, position, index of the variant and centimorgan.
* each additional column specifies the probability of copying a target haplotype from a group of individuals.

* Ex: a column named ```NA20900_G1_0``` stores the probability of copying the haplotype ```0``` of the individuals ```NA20900``` from the group ```G1```.
* Ex: a column named ```NA20900_HOLE_1``` stores the probability of copying the haplotype ```1``` of the individuals ```NA20900``` from the group ```H``` (unrelated individuals).



</br>
  &nbsp


#
## Installation

A static version of THORIN ready to run is provided here: ```./bin/thorin_static```



To build THORIN from source, follow the steps below.

### System requirements.
THORIN is a C++ tool. In order to compile, we require a modern Linux operating system and a version of GCC > 4.4. We recommend to use the latest available version for your system.

For example running the following instruction on Ubuntu 20.04 focal:

	sudo apt install build-essential
	
will install the GNU g++ compiler version 9.2. To check the version of your g++ compiler, simply run:

	g++ --version
	
### Required libraries.
THORIN requires several libraries installed on the system. Here we assume most of the libraries are available, and we focus on two main libraries:

- HTSlib version >= 1.7: A C library for reading/writing high-throughput sequencing data.
- BOOST version >= 1.65: A set of peer-reviewed portable C++ source libraries.


#### HTSlib

Building HTSlib is straightforward and does not require root privileges. Please refer to the HTSlib documentation for complete details. Here we provide a basic script to install HTSlib v1.16:

	wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2
	tar -xf htslib-1.16.tar.bz2
	mv htslib-1.16 htslib
	cd htslib
	autoheader; autoconf; ./configure; #optional
	make

After this, you’ll find the libhts.a inside your current directory and the include files inside subdirectory: ./include/


#### Boost

As THORIN only requires few of the boost libraries, we’ll build the smallest possible boost build, without requiring root privileges. Please refer to the Boost installation instructions for complete details. Here we provide a basic script to the minimal build of Boost v1.73.0 required to run THORIN:

	wget https://boostorg.jfrog.io/artifactory/main/release/1.73.0/source/boost_1_73_0.tar.bz2
	tar --bzip2 -xf boost_1_73_0.tar.bz2
	cd boost_1_73_0
	./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost #where ../boost is your custom boost installation prefix
	./b2 install
	cd ../boost #change this to the folder you used as --prefix for the bootstrap script

After this, you will also find a copy of the Boost headers in the include/ subdirectory of the installation prefix (our current directory). The Boost static libraries (libboost_iostreams.a, libboost_program_options.a and libboost_serialization.a) are found in the subfolder ./lib/ of your installation prefix.


3. Additional libraries

Make sure that the following standard library flags can be used by g++ on your system:

    -lz,-lbz2 and -llzma for reading/writing compressed files.
    -lm for basic math operations.
    -lpthread for multi-threading

You can do so by checking the outcome of the following commands:

	locate -b '\libz.so'
	locate -b '\libbz2.so'
	locate -b '\liblzma.so'
	locate -b '\libm.so'
	locate -b '\libpthread.so'
	locate -b '\libcurl.so'







</br>
  &nbsp
  
 #
## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

