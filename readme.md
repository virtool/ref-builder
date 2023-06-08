# Reference Builder

Automatically build and maintain reference sets of pathogen genome sequences.

## TODO

* Do we allow 

## Intro

### Background

Online genomic sequence databases are important tools for the dissemination of genomic information. This information is especially for developing and performing diagnostic assay based on pathogen genome sequences.

GenBank is the most well known repository of genomic sequence data. It accepts a broad array of sequence data and is an invaluable resource for research. However, it is not strongly curated, structured, or validated. This makes it flexible, but allows for duplications and incorrect or incomplete information. For these reasons, Genbank is not an ideal resource for making diagnostic decisions.

### Mission

This framework provides a way to build and maintain reference sequence datasets using GitHub. This will support collaboration on transparent, shared reference sets for genomics-based diagnostics.

## Features

### Automation

* Genbank is regularly scanned for new virus isolates. These are submitted as PRs with human-readable comments and checks that can be interpreted by curators.
* Taxonomy DB is regularly scanned regularly for taxonomy IDs for viruses that do not already have one.

## Guide

### Submit a Virus

Submit an issue with the following information:

* Virus Name
* Isolate Name

[Submit a Virus](https://github.com/virtool/ref-builder/issues)

Submissions are started as GitHub issues. Users provide the taxonomy isolates names, 

### Submit an Isolate

Submit an issue with the following information:

* Isolate Type (eg.  Clone,  Extract, Isolate, Strain, Variant)
* Isolate Name (eg. A, BC-1, 2022.10932)

[Submit an Isolate](https://github.com/virtool/ref-builder/issues)

## FAQ

**Can I used Ref Builder with GitLab, BitBucket or others?**

No. We rely on GitHub's features to make Ref Builder work. Please contact us if you'd
like to contribute to support for other platforms.

**Can I submit a new virus with no taxonomy ID?**

Yes. Ref Builder will scan Taxonomy DB for a taxonomy ID for your virus. It will create an issue after a while if your virus still
doesn't have a taxonomy ID.
