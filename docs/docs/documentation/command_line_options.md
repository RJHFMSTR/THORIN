---
layout: default
title: Command line options
nav_order: 2
parent: Documentation
---

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-\-seed             | INT     | 15052011 | Seed of the random number generator  |
| \-T \[ \-\-thread \] | INT     | 1        | Number of thread used|

#### Input files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Genotypes of focal and reference individuals in VCF/BCF/ format |
| \-H \[\-\-hole \]| STRING  | NA       | Genotypes of unrelated individuals in VCF/BCF/ format  |
| \-G \[\-\-group \]| STRING  | NA       | Group file listing focal and reference individuals |
| \-M \[\-\-map \]     | STRING  | NA       | Genetic map  |
| \-R \[\-\-region \]  | STRING  | NA       | Target region  |


#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Sum of HMM IBD probabilities per group TXT/VCF/BCF format |
| \-\-ibd     | STRING  | NA       | IBD segments per group  |
| \-\-phasing  | STRING  | NA       | IBD-based scaffold in VCF/BCF format |
| \-\-log              | STRING  | NA       | Log file  |



