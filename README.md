# **SARS-CoV-2 Sequencing Resources**

This document repository is meant to serve as the start of a crowd-sourced collection of information, documentation, protocols and other resources for public health laboratories intending to sequence SARS-CoV-2 coronavirus samples in the coming weeks. This is admittedly a limited first draft, but will continued to collate useful information as additional protocols, tools, and resources are added, and as best practices are identified. While some of the resources here are directed specifically to US state and local public health laboratories in support of diagnostic testing, sequencing and response, we hope that this is a useful resource for the global laboratory community, as we respond to this pandemic threat.

This collection is maintained and curated by [Duncan MacCannell](https://www.twitter.com/dmaccannell) from the Office of Advanced Molecular Detection ([AMD](https://www.cdc.gov/amd)) at the Centers for Disease Control and Prevention ([CDC](https://www.cdc.gov)). Please feel free to suggest additions, edits, clarifications and corrections -- either by posting an issue, filing a pull request or by contacting me directly by email or twitter.  In the meantime, I'll continue to add and mirror useful resources here as they become available.

*Disclaimer*
> The findings and conclusions in this document and the attendant repository are those of the author and do not necessarily represent the official position of the Centers for Disease Control and Prevention. Use of trade names is for identiﬁcation only and does not imply endorsement by the Centers for Disease Control and Prevention or by the U.S. Department of Health and Human Services.


---


# **Sequencing Protocols**

## _**1. OXFORD NANOPORE**_
The following sequencing protocols, checklists and job-aids are primarily designed for the Oxford Nanopore [MinION](https://nanoporetech.com/products/minion), and have been kindly shared by research groups throughout the world (please see individual protocols for attribution and citing purposes). Even so, most of these protocols should scale to larger ONT instruments without significant modifications.

### a) CDC NCIRD/DVD ONT Sequencing Protocol
This protocol was developed, tuned and validated by the Viral Discovery laboratory at CDC/NCIRD, where it was used to generate the first 16 SARS-CoV-2 genome sequences from the United States. In practice, it has been used for situations with a relatively low or predictable volume of samples, and is often used in conjunction with Sanger-based tiling to resolve any potential sequencing or assembly issues.
- Release pending CDC clearance.

### b) [ARTIC Network nCoV-2019 Sequencing Protocol](https://artic.network/ncov-2019)
This protocol was developed and released by the fine folks at [ARTIC Network](https://artic.network), and was subsequently refined based on comments from [Itokawa et al](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v1.full.pdf), which identified potential issues and proposed an alternate L18 primer.

- [Sequencing protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w) / [Single sample sequencing protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-single-sample-bdbfi2jn)

- [Stepwise simplified protocol from ONT*](./protocols/ONT-COVID-19_Tiling)

- Primer schemes: [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1) / [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) [(ref)](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v1.full.pdf)

- Integrated bioinformatics (RAMPART) - documentation below under bioinformatics methods.

### c) Doherty Institute VIDRL Sequencing Protocols
The Victorian Infectious Diseases Reference Laboratory ([VIDRL](https://www.vidrl.org.au/)) at the [Peter Doherty Institute for Infection and Immunity](https://www.doherty.edu.au/) released two protocols for the ONT MinION, which they successfully used to sequence early Australian SARS-CoV-2 samples.
- [Native RNA Sequencing Protocol](./protocols/ONT-Native_RNA) [(ref)](https://www.biorxiv.org/content/10.1101/2020.03.05.976167v1.full.pdf)

- [SISPA Sequencing Protocol](./protocols/ONT-SISPA)

## _**2. ILLUMINA**_
### a) Illumina Nextera Flex Enrichment Sequencing Protocol
Illumina's Research and Development group has recently developed and validated a custom, research use only (RUO) enrichment sequencing strategy based on their Nextera Flex chemistry.
- Release imminent.

### b) SARS-CoV-2 Enrichment Sequencing by Spiked Primer MSSPE (NIAID/UCSF/CZBioHub)
The NIAID laboratory team in Cambodia, in collaboration with UCSF, CZBioHub and IPC, has released a metagenomic sequencing with spiked primer enrichment (MSSPE) protocol for SARS-CoV-2. The protocol is available on [protocols.io](https://protocols.io).
- [MSSPE Protocol](https://www.protocols.io/view/sars-cov-2-enrichment-sequencing-by-spiked-primer-bc36iyre)

### c) Illumina Shotgun Metagenomics Sequencing Protocol
Illumina's technical note on sequencing coronavirus samples using a comprehensive metagenomic sequencing approach was one of the earlier protocols released for SARS-CoV-2, and remains an effective option for shotgun sequencing.
- [Application Note](https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/ngs-coronavirus-app-note-1270-2020-001.pdf)

## _**3. TANDEM**_
### a) SARS-CoV-2 Parallel Sequencing by Illumina and ONT (UWMadison ZEST)
Staff and students from Thomas Friedrich and Dave O'Connor's laboratories at UWMadison have put together a tandem sequencing protocol and bioinformatic workflow that incorporates Illumina and ONT sequence. While this may be overkill for routine or high-throughput public health purposes, the necessary protocols, scripts and documentation are available here.
- [Protocol](https://openresearch.labkey.com/wiki/ZEST/Ncov/page.view?name=SARS0CoV-2%20Deep%20Sequencing)
- [Github resources](https://github.com/katarinabraun/SARS-CoV-2_sequencing)

## _**4. SANGER**_
### a) CDC NCIRD/DVD Sanger Tiling
An elegant approach from a more civilized age. The Viral Discovery laboratory team in CDC/NCIRD/DVD has used conventional Sanger sequencing to refine and complement betacoronavirus sequencing on next generation platforms.
- Release pending CDC clearance.
---


# **Bioinformatics Tools, Scripts and Workflows**

### 1. CDC NCIRD/DVD Bioinformatics SOPs
This section describes the basic bioinformatic workflow that the Viral Discovery laboratory in NCIRD, and other teams at CDC use for quality assessment, assembly and comparison of coronavirus sequences. IRMA, the Iterative Refinement Meta-Assembler developed by CDC's Influenza Division for routine influenza surveillance, has recently been updated to support both ebolavirus and coronavirus assembly tasks. While IRMA isn't used for all SARS-CoV-2 assemblies at CDC, it is a powerful tool for complex or problematic samples and datasets.
- Release of scripts and other tools is pending CDC clearance.
- [IRMA: Iterative Refinement Meta-Assembler](https://wonder.cdc.gov/amd/flu/irma) is available here.

### 2. CLCbio Genomics Workbench
QIAGEN CLCbio will be releasing a set of CLCbio SARS-CoV-2 workflows that have been optimized for Illumina, ONT and IonTorrent sequencing platforms.  More information, download instructions and tutorials will be available within the next several days.
- SARSCoV2 Assembly/Analysis Tutorial and FAQs (Illumina/ONT/Ion data) - TBD

### 3. ARTIC Network Bioinformatics
The [ARTIC Network](https://artic.network) has released detailed instructions on how to setup and configure the conda environment needed to run their analysis pipelines. These are complete bioinformatic workflows, including runtime visualization, basecalling, mapping/assembly and reporting in a single, portable environment. The [artic-nCoV2019](https://github.com/artic-network/artic-ncov2019) repo includes source code and build instructions for a custom RAMPART configuration.  Additional instructions and documentation are available below.

- [Bioinformatics Environment Setup/Configuration](https://artic.network/ncov-2019/ncov2019-it-setup.html)

- [SARS-CoV-2 Bioinformatics SOP](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)

- [RAMPART Runtime Documentation](https://artic.network/ncov-2019/ncov2019-using-rampart.html)


---

# **Quality Management**
This section will describe best practices for laboratory and bioinformatic quality assurance, including preflight checks for sequence and metadata submission to public repositories.
- **to do**: I'd love for people to help describe their actual QC processes.

---


# **Submitting to Public Sequence Repositories**

### Sequence naming conventions for public repositories
We are proposing simplified naming conventions for sequences submitted to GISAID and NCBI from US public health and clinical laboratories.

**COUNTRY** | /  **STATE-LAB-SAMPLE**  / | **YEAR**
--- | --- | ---
`USA` | ``/ CA-CDPH-S1 /`` | `2020`
`USA` | ``/ UT-01 /`` | `2020`
`USA` | ``/ AZ-1045 /`` | `2020`
`USA` | ``/ WA-UW-316 /`` | `2020`

> The proposed convention is as follows: 1) country (USA). 2) The middle sample identification cell should include two-letter state (eg: CA), an abbreviated identifier for the submitting lab (eg: CDPH), as desired, and a unique sequence identifier (eg:01, S01, 454, ...), with all three terms separated by hyphens.




For states with only one submitting laboratory (which should be most), the identifier for the submitting laboratory may be omitted, resulting in a simple, state-level identifier such as `USA/UT-573/2020`.

** These recommendations are roughly compatible with existing submissions to GISAID and NCBI, but are completely open for debate.*

### Recommended formatting and criteria for sample metadata

#### [NCBI SARS-CoV-2 Sequences](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/) [(submit)](https://www.ncbi.nlm.nih.gov/WebSub/)
The National Center for Biotechnology has established a custom landing page for SARS-CoV-2 sequences and data, and is working to develop streamlined submission processes for Genbank and SRA.  For the time being, we suggest basing metadata and submission formatting on GISAID EpiCoV, which tends to be more comprehensive and structured. We will develop specific guidance for NCBI submissions.

#### [GISAID EpiCoV](https://www.gisaid.org/)

The GISAID EpiCoV Public Access repository is based on existing submission processes and data structures for large-scale influenza surveillance (GISAID EpiFlu). As such, submitters to EpiCoV will discover that several of the required metadata submission fields may be problematic.  Nonetheless, a number of laboratories have been submitting sequences with the following:

METADATA FIELDS (GISAID) | GUIDANCE
--- | ---
`Virus name`| USA/FL-103/2020 (see above)
`Accession ID` |
`Type` | betacoronavirus
`Collection date`| YYYY-MM-DD
`Location` | USA / State / County?
`Additional location information` |
`Host` | Human
`Additional host information` |
`Gender` | (no guidance)
`Patient age`| (no guidance, could be binned)
`Patient status` | (no guidance)
`Specimen source` | (free text)
`Outbreak detail` | omit
`Last vaccinated` | omit
`Treatment` | omit

At a minimum, we suggest that samples be submitted with `collection date` `location` `host` information attached.  `location`, `host`, `gender` `patient age` are all required fields, and several of them likely constitute personally-identifiable information. While they cannot be left blank for submission, you can submit the record successfully (in both single or batch mode) by entering "**unknown**".

Note that for GISAID submissions, users must register for an account, and must successfully submit a single submission before being granted access to the bulk submission template and interface.

---


# **Linking Sequence Accessions**
For data linkage, we are proposing the following template, as a simple, lightweight line list of tab-separated values. If this consensus recommendation for data linkage is acceptable, a preformatted .TSV will be made available. We recognize that not all samples sent for sequencing have a PUID associated.

STATE | GISAID ID | GENBANK ID | COLLECTION DATE | PUID/COVID-ID
--- | --- | --- | --- | --- 

In this simple proposed schema, GISAID ID or GENBANK ID and COLLECTION DATE are required fields, and our hope is to maximize PUID completion.  All accession numbers, including PUID should be entered without any superfluous text or annotation.

---


# **Other Useful References and Resources**

### [Nextstrain hCoV-19](https://www.nextstrain.org/ncov)
### [Virological Novel 2019 Coronavirus forum](http://virological.org/c/novel-2019-coronavirus)
### [CNCB 2019 Novel Coronavirus Resource (2019nCoVR)](https://bigd.big.ac.cn/ncov?lang=en)



---

Updated: 20200316 [@dmaccannell](https://www.twitter.com/dmaccannell)
