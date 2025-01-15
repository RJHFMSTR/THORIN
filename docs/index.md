---
layout: default
title: Home
nav_order: 1
description: "THORIN is a tool for mapping Identity-By-Descent (IBD) and inferring the Parent-of-Origin of alleles."
permalink: /

---

![](assets/images/logo_name.png?raw=true)

<!---
# THORIN
{: .fs-9 .fw-500 }
-->

<!---
**T**arget **H**aplotype **OR**igin **IN**ference version **1.2**
{: .fs-5 }
-->

---

## About

THORIN is a tool for mapping Identify-By-Descent (IBD) and inferring the Parent-of-Origin of alleles.

## News

{: .new }
> **Version `1.2.0` is now available!**
<!--- > See [the CHANGELOG](https://github.com/odelaneau/shapeit5/blob/main/docs/CHANGELOG.md) for details.
-->

## Citation

If you use THORIN in your research, please cite the following paper:

[Hofmeister, R.J. et al. Parent-of-Origin inference for biobanks. Nature Communication 2022. https://doi.org/10.1038/s41467-022-34383-6](https://www.nature.com/articles/s41467-022-34383-6)

[Hofmeister, R.J. et al. Parent-of-Origin inference and its role in the genetic architecture of complex traits: evidence from ∼220,000 individuals](https://www.medrxiv.org/content/10.1101/2024.12.03.24318392v1)

---

[Get started now](#getting-started){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 .mx-auto }
[View source code on GitHub](https://github.com/rjhfmstr/thorin){: .btn .fs-5 .mb-4 .mb-md-0 }



---

## Getting started

- [See documentation](https://rjhfmstr.github.io/THORIN/docs/documentation)

---

## Description

THORIN is a tool to map Identity-by-Descent between a focal individual and reference individuals that is particularly well suited for inferring the parent-of-origin of alleles, by using surrogate parent groups as reference individuals. It is composed of the following options:

- **IBD probability per variant site**. THORIN can output the probability of sharing IBD at each variant sites.
- **IBD probability per variant site, VCF/BCF format**. The output file can be formatted as indexable VCF/BCF file to allow easy extraction of individuals or genomic region of interests.
- **IBD segments**. THORIN can output shared IBD haplotype segments.
- **IBD scaffold**. THORIN can provide an IBD-based scaffold for re-estimating haplotype, both for correcting long-range intra-chromosomal phasing errors and for performing inter-chromosomal phasing.


[IBD per variant site](https://rjhfmstr.github.io/THORIN/docs/documentation/inputs_and_outputs.html#1-ibd-per-variant-site){: .btn .btn-blue }
[IBD probability per variant site, VCF/BCF format](https://rjhfmstr.github.io/THORIN/docs/documentation/inputs_and_outputs.html#2-ibd-per-variant-site-variant-call-format){: .btn .btn-blue }
[IBD segments](https://rjhfmstr.github.io/THORIN/docs/documentation/inputs_and_outputs.html#3-ibd-segments){: .btn .btn-blue }
[IBD scaffold](https://rjhfmstr.github.io/THORIN/docs/documentation/inputs_and_outputs.html#4-ibd-scaffold){: .btn .btn-blue }

---

## About the project

THORIN is developed by Robin Hofmeister, Theo Cavinato and Olivier Delaneau.

### License

THORIN is distributed with [MIT license](https://github.com/RJHFMSTR/THORIN/blob/main/LICENSE).

### Organisations

<div class="d-flex justify-content-around">
  <div class="p-5"><a href="https://www.unil.ch/index.html"><img src="assets/images/lausanne_logo.jpg" align="right" alt="unil" style="height:50px"></a></div>
  <div class="p-5"><a href="https://www.sib.swiss/"><img src="assets/images/sib_logo.jpg" align="right" alt="sib" style="height:50px"></a></div>
  <div class="p-5"><a href="https://www.snf.ch/en/Pages/default.aspx"><img src="assets/images/snf.gif" align="right" alt="snf" style="height:50px"></a></div>
</div>

### Contributing

THORIN is an open source project and we very much welcome new contributors. When contributing to our repository, please first discuss the change you wish to make via issue,
email, or any other method with the owners of this repository before making a change.
#### Thank you to the contributors of THORIN!

{% if site.github.contributors %}
  {% for contributor in site.github.contributors %}
    <a href="{{ contributor.html_url }}">
      <img src="{{ contributor.avatar_url }}" width="32" height="32" alt="{{ contributor.login }}" />
    </a>
  {% endfor %}
{% else %}
  <p>No contributors found.</p>
{% endif %}


---

We thank the [Just the Docs](https://github.com/just-the-docs/just-the-docs) developers, who made this awesome theme for Jekyll.



