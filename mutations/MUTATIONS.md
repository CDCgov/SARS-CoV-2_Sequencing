## Common disruptive mutations observed in SARS-CoV-2 genome sequences

A list of mutations observed by public health or academic institutions engaged in sequencing SARS-CoV-2 clinical specimens. Listed mutations are submitted by the community and come with no guarantees or endorsements. Mutations are provided in both tabular and VCF formats.  

See also https://github.com/W-L/ProblematicSites_SARS-CoV2

### Intended use

1. Facilitate sample review by public health laboratories.
1. Aid review of assemblies flagged by QC processes due to presumed disruptive mutations.
1. Compile problematic sites that may impact downstream bioinformatic analyses.
1. Xxxxxx

### Inclusion criteria

Before adding, submitters are expected to confirm the validity of mutations to ensure high confidence according to the following suggested criteria:  
1. Read coverage depth >= Xx
1. Xxxxxx
1. Xxxxxx


### How to contribute mutations

1. Fork repo  
1. Add entries to 'mutations.tsv'
1. Update VCF file with 'some-script.sh'
1. Commit changes and submit a pull request

### Table column definitions

| Header         | Description                    |
|----------------|--------------------------------|
|Chrom		|Name of the reference sequence (Wuhan-1)|
|Gene			|Annotated gene in the reference|
|Pos 			|Base position in the reference|
|Mutation	|Mutation change|
|Description |Human-readable description of mutation|
|Example	|GenBank and/or SRA accession of example genome with confirmed mutation|
|Library	|Primer set or library method information|
|Platform	|Sequencing platform information|
|Analysis	|Bioinformatic pipeline information (version)|
|Lineage	|Affected PANGO lineages (Pangolin version)|
|Submitter	|Name and affiliation of mutation submitter|
|Date	|Date added to the table|

### VCF column definitions (v4.? specification)

| Header         | Description                    |
|----------------|--------------------------------|
|CHROM           | Name of the reference sequence |
|POS             | 1-based position of the variation on the reference |
|ID              | NA? |
|REF             | Reference base |
|ALT             | List of alternative alleles at the position (IUPAC ambiguity code) |
|QUAL            | NA? |
|FILTER          | ???? |
|INFO            | ???? |
