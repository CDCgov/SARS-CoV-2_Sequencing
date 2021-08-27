## Common disruptive mutations observed in SARS-CoV-2 genome sequences

A list of mutations observed by public health or academic institutions engaged in sequencing SARS-CoV-2 clinical specimens. Listed mutations are submitted by the community and come with no guarantees or endorsements. Mutations relative to reference [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) are provided in both table and VCF formats as described below.  

See also: https://github.com/W-L/ProblematicSites_SARS-CoV2

### Intended use

1. Facilitate sample review by public health laboratories.
1. Aid review of assemblies flagged by QC processes due to presumed disruptive mutations.
1. Compile problematic sites that may impact downstream bioinformatic analyses.
1. Xxxxxx

### Inclusion criteria

Before adding, submitters are expected to confirm the validity of mutations to ensure high confidence according to the following suggested criteria for each platform:  

<details>
 <summary>Illumina MiSeq/NextSeq</summary>

1. Read coverage depth >= Xx  
1. Xxxxxx
</details>

<details>
 <summary>Oxford Nanopore (and ClearLabs)</summary>

1. Read coverage depth >= Xx
1. Xxxxxx
</details>

<details>
 <summary>PacBio</summary>

1. Read coverage depth >= Xx
1. Xxxxxx
</details>


### How to contribute mutations

1. Fork this repository ([Do you need a GitHub account?](https://docs.github.com/en/get-started/signing-up-for-github/signing-up-for-a-new-github-account))  
1. Copy `./import/mut_template.csv` and name your file with an unique identifier using format:
    [YYYYMMDD].[YourName].[0-9]{3}.tsv that can be used to link entries between the table and changelog (eg. `./20210812.Weigand.123.csv`).    
1. Enter mutation information according to the column descriptions below, entering 'NA' for any missing attributes except those marked as required. You may add multiple entries to this file.
1. Save your new mutation file to `./import`.
1. Commit your changes.
1. Submit a pull request.
1. Updated versions of `SCV2_mutations_latest.csv` and `SCV2_mutations_latest.vcf` will be generated automatically when your pull request is approved.

### Table column definitions

| Header         | Description                    |
|----------------|--------------------------------|
|Chrom *		|Name of the reference sequence [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2)|
|Gene	*		|Annotated gene in the reference|
|Position *			|Base position in the reference|
|Reference *	|Reference sequence |
|Mutation * |Mutation sequence change (use '-' for deletion)|
|Description * |Human-readable description of mutation|
|Example	|GenBank and/or SRA accession of example genome with confirmed mutation|
|Library	|Primer set or library method information|
|Platform	|Sequencing platform information|
|Analysis	|Bioinformatic pipeline information (version)|
|Lineage	|Affected PANGO lineages (Pangolin version)|
|Submitter *	|Name and affiliation of mutation submitter|
|Date *	|Date added to table|

(\* required)

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
