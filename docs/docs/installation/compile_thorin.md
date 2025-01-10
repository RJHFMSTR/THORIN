---
layout: default
title: Compile THORIN
parent: Build from source
grand_parent: Installation
permalink: /docs/installation/build_from_source/compile_thorin
---
# Compile SHAPEIT5
{: .no_toc .text-center }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

## Compile THORIN
Download the last version of the THORIN code using:
<div class="code-example" markdown="1">
```bash
git clone https://github.com/RJHFMSTR/THORIN.git
```
</div>


Navigate to the downloaded folder using `cd THORIN`.

You'll find there a folder containing the software with the structure:

- `bin`: folder for the compiled binary.
- `obj`: folder with all binary objects.
- `src`: folder with source code.
- `makefile`: Makefile to compile the program.

In order to compile the THORIN tool, you need to edit the makefile at lines so that the following variables are correctly set up (look at the paths already there for an example):

- `HTSSRC`: path to the root of the HTSlib library, the prefix for HTSLIB_INC and HTSLIB_LIB paths.
- `HTSLIB_INC`: path to the HTSlib header files.
- `HTSLIB_LIB`: path to the static HTSlib library (file `libhts.a`).
- `BOOST_INC`: path to the BOOST header files (usually `/usr/include`). 
- `BOOST_LIB_IO`: path to the static BOOST iostreams library (file `libboost_iostreams.a`). 
- `BOOST_LIB_PO`: path to the static BOOST program_options library (file `libboost_program_options.a`). 

If installed at the system level, static libraries (*.a files) can be located with this command:

<div class="code-example" markdown="1">
```bash
locate libboost_program_options.a libboost_iostreams.a libhts.a
```
</div>

Once all paths correctly set up, proceed with the compilation using `make`. The binary can be found in the `bin/` folder.


