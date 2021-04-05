# NCBI TaxID Counter
This project is part of a collection of codes developed for my Undergraduate Thesis.

## Dependencies
* Biopython - Bio.Entrez package
* Matplotlib - venn2, venn2_circle, pyplot
* Modules - csv, os and sys
* python 3.6+

## Features
* It takes the result files from Metagenome Taxonomy Classification Tools and sorts the TaxIDs between the taxonomy ranks (Superkingdom, Kingdom, Phylum, Class, Order, Family, Genus and Species) and the groups used by the programs (Bacteria, Archaea, Viruses, Fungi, Protozoa)

## Usage
### Arguments
It can handle as many files as need, each being processed individually.
#### Example:
```bash
./taxid_counter.py (-compare) krona_file metacache_file centrifuge_file [...] 
```
### Output
- The output is a CSV file where the rows represent same rank and the columns same group. The name is the original file name + _sorted.csv" 
- It also outputs Venn diagrams in SVG format for Bacteria, Archaea, Viruses, Fungi and Protozoa with total TaxIds and only species.
