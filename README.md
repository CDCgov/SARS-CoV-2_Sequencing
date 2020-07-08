# **SARS-CoV-2 Sequencing Resources**

This document repository is meant to serve as the start of a crowd-sourced collection of information, documentation, protocols and other resources for public health laboratories intending to sequence SARS-CoV-2 coronavirus samples in the coming weeks. This is admittedly a limited first draft, but will continued to collate useful information as additional protocols, tools, and resources are added, and as best practices are identified. While some of the resources here are directed specifically to US state and local public health laboratories in support of diagnostic testing, sequencing and response, we hope that this is a useful resource for the global laboratory community, as we respond to this pandemic threat.

This collection is maintained and curated by [Duncan MacCannell](https://www.twitter.com/dmaccannell) from the Office of Advanced Molecular Detection ([AMD](https://www.cdc.gov/amd)) at the Centers for Disease Control and Prevention ([CDC](https://www.cdc.gov)). Please feel free to suggest additions, edits, clarifications and corrections -- either by posting an issue, filing a pull request or by contacting me directly by email or twitter.  In the meantime, I'll continue to add and mirror useful resources here as they become available.

## INDEX
- [Sequencing Protocols](#sequencing)
- [Bioinformatic Tools, Scripts and Workflows](#bioinformatics)
- [Quality Management](#quality)
- [Submitting to Public Sequence Repositories](#submitting)
- [Linking Sequence Accessions](#linking)
- [Other Useful References and Resources](#resources)
- [Notices and Disclaimers](#notices)


*Disclaimer*
> The findings and conclusions in this document and the attendant repository are those of the author and do not necessarily represent the official position of the Centers for Disease Control and Prevention. Use of trade names is for identiﬁcation only and does not imply endorsement by the Centers for Disease Control and Prevention or by the U.S. Department of Health and Human Services.


---


# **Sequencing Protocols**<a id="sequencing" />

## _**1. OXFORD NANOPORE**_
The following sequencing protocols, checklists and job-aids are primarily designed for the Oxford Nanopore [MinION](https://nanoporetech.com/products/minion), and have been kindly shared by research groups throughout the world (please see individual protocols for attribution and citing purposes). Even so, most of these protocols should scale to larger ONT instruments without significant modifications.

### a) CDC NCIRD/DVD ONT Sequencing Protocol
This protocol was developed, tuned and validated by the Viral Discovery laboratory at CDC/NCIRD, where it was used to generate the first 16 SARS-CoV-2 genome sequences from the United States. In practice, it has been used for situations with a relatively low or predictable volume of samples, and is often used in conjunction with Sanger-based tiling to resolve any potential sequencing or assembly issues.
- [Singleplex and Multiplex Protocols for ONT/Illumina](./protocols/CDC-Comprehensive)

### b) [ARTIC Network nCoV-2019 Sequencing Protocol](https://artic.network/ncov-2019)
This protocol was developed and released by the fine folks at [ARTIC Network](https://artic.network), and was subsequently refined based on comments from [Itokawa et al](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v1.full.pdf), which identified potential issues and proposed an alternate L18 primer.

- [Sequencing protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-bbmuik6w) / [Single sample sequencing protocol](https://www.protocols.io/view/ncov-2019-sequencing-protocol-single-sample-bdbfi2jn)

- [Stepwise simplified protocol from ONT*](./protocols/ONT-COVID-19_Tiling)

- Primer schemes: [V1](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V1) / [V2](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V2) [(ref)](https://www.biorxiv.org/content/10.1101/2020.03.10.985150v1.full.pdf) / [V3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3)

- [ARTIC on Illumina - Complete Walk-through - Grubaugh/Andersen/Loman labs](https://docs.google.com/document/d/1PilT4w5jHO-ROsE8TL5WBGa0wSCdTHAsNl1LIOYiTgk/mobilebasic)

- Integrated bioinformatics (RAMPART) - documentation below under bioinformatics methods.

### c) Doherty Institute VIDRL Sequencing Protocols
The Victorian Infectious Diseases Reference Laboratory ([VIDRL](https://www.vidrl.org.au/)) at the [Peter Doherty Institute for Infection and Immunity](https://www.doherty.edu.au/) released two protocols for the ONT MinION, which they successfully used to sequence early Australian SARS-CoV-2 samples.
- [Native RNA Sequencing Protocol](./protocols/ONT-Native_RNA) [(ref)](https://www.biorxiv.org/content/10.1101/2020.03.05.976167v1.full.pdf)

- [SISPA Sequencing Protocol](./protocols/ONT-SISPA)

## _**2. ILLUMINA**_
### a) Illumina Nextera Flex Enrichment Sequencing Protocol
Illumina's Research and Development group has recently developed and validated a custom, research use only (RUO) enrichment sequencing strategy based on their Nextera Flex chemistry.
- [Application Note](./protocols/ILMN-Nextera_Flex_Enrichment/ngs-enrichment-coronavirus-app-note-FINAL.pdf)

### b) SARS-CoV-2 Enrichment Sequencing by Spiked Primer MSSPE (NIAID/UCSF/CZBioHub)
The NIAID laboratory team in Cambodia, in collaboration with UCSF, CZBioHub and IPC, has released a metagenomic sequencing with spiked primer enrichment (MSSPE) protocol for SARS-CoV-2. The protocol is available on [protocols.io](https://protocols.io).
- [MSSPE Protocol](https://www.protocols.io/view/sars-cov-2-enrichment-sequencing-by-spiked-primer-bc36iyre)

### c) Illumina Shotgun Metagenomics Sequencing Protocol
Illumina's technical note on sequencing coronavirus samples using a comprehensive metagenomic sequencing approach was one of the earlier protocols released for SARS-CoV-2, and remains an effective option for shotgun sequencing.
- [Application Note](https://www.illumina.com/content/dam/illumina-marketing/documents/products/appnotes/ngs-coronavirus-app-note-1270-2020-001.pdf)

### d) SARS-CoV-2 and related virus sequencing with capture enrichment (Broad Institute)
The Sabeti lab, at the Broad Institute, released a probe set for comprehensive whole-genome capture of SARS-CoV-2 and respiratory-related viruses (human-infecting coronaviruses, HRSV, HMPV, HPIVs, Human mastadenovirus A-G, Enterovirus A-E, Rhinovirus A/B/C, influenza A/B/C).
The probe set is available as V-Respiratory on the [probe designs page](https://github.com/broadinstitute/catch/tree/master/probe-designs) of the CATCH repository.
It was initially released in January, 2020 and most recently updated in March, 2020.
Probes can be ordered from Twist Bioscience; we have used the protocol for [Twist custom panels](https://www.twistbioscience.com/products/ngs#product-featured-2911) with slight modifications for low input Nextera XT libraries.
- [Probe set](https://github.com/broadinstitute/catch/tree/master/probe-designs)
- [Design approach and prior results](https://www.nature.com/articles/s41587-018-0006-x)

### e) ARTIC on Illumina
A number of different laboratories have implemented derivatives of the ARTIC amplicon scheme on Illumina.
- [Protocol from the Grubaugh lab at Yale (gdoc)](https://docs.google.com/document/d/1PilT4w5jHO-ROsE8TL5WBGa0wSCdTHAsNl1LIOYiTgk/edit)
- [Virgina DCLS has a detailed ARTIC-on-Illumina protocol, including worksheets, flowcharts and helpful hints.](./protocols/ILMN-VA_DCLS-ARTIC_Illumina)
- Joel Sevinsky and the intrepid team at [StaPH-B](https://www.staphb.org) have adapted ARTICv3 for Illumina for public health laboratories that are already configured for PulseNet sequencing.
      - [SARS-CoV-2 Sequencing on Illumina MiSeq using ARTIC Protocol: Part 1 - Tiling PCR](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bfefjjbn)
      - [SARS-CoV-2 Sequencing on Illumina MiSeq using ARTIC Protocol: Part 2 - Illumina DNA Flex](https://www.protocols.io/view/sars-cov-2-sequencing-on-illumina-miseq-using-arti-bffyjjpw) 

#### Modifications, Improvements and Derivatives of ARTIC
- [Freed and Silander developed a modified ARTIC protocol with 1200bp inserts](https://www.protocols.io/view/ncov-2019-sequencing-protocol-rapid-barcoding-1200-bfwnjpde)([primers](https://docs.google.com/spreadsheets/d/1M5I_C56ZC8_2Ycgm9EFieVlVNqxsP7dXAnGoBZy3nDo/edit#gid=755704891))
- Gohl and collagues at UMN describe an elegant tailed amplicon approach that enables ARTIC without the need for ligation or tagmentation steps ([biorXiv](https://www.biorxiv.org/content/10.1101/2020.05.11.088724v1))([protocols.io](https://www.protocols.io/view/sars-cov-2-tailed-amplicon-illumina-sequencing-bge5jtg6))
- Sorensen, Karst and Knutsson from Alborg University have described a tailed long amplicon derivative of ARTIC ([protocols.io](https://www.protocols.io/view/aau-ncov-2019-tailed-long-amplicon-sequncing-bfc3jiyn))

## _**3. PACIFIC BIOSCIENCES**_
PacBio maintains a [COVID-19 landing page](https://www.pacb.com/research-focus/microbiology/COVID-19-sequencing-tools-and-resources/) with updated resources on existing protocols and [SARS-CoV2-2 assay development recommendations](https://www.pacb.com/wp-content/uploads/Customer-Collaboration-SARS-CoV-2-Virus-Sample-Preparation-Suggestions-for-PacBio-SMRT-Sequencing.pdf).

### Barcoding options
* [Barcoding with Barcoded Overhang Adapters](https://www.pacb.com/wp-content/uploads/Procedure-Checklist-Preparing-SMRTbell-Libraries-using-PacBio-Barcoded-Overhang-Adapters-for-Multiplexing-Amplicons.pdf)
* [Barcoding with M13 primers tailed with PacBio Barcodes](https://www.pacb.com/wp-content/uploads/Procedure-Checklist-Preparing-SMRTbell-Libraries-using-PacBio-Barcoded-M13-Primers-for-Multiplex-SMRT-Sequencing.pdf)

a) [Sinai 1.5kb and 2kb tiled amplicon protocol](https://www.medrxiv.org/content/10.1101/2020.04.08.20056929v1.full.pdf) and [PacBio barcoding options](https://www.pacb.com/wp-content/uploads/Sinai-SARS-CoV-2-Primers-with-PacBio-Barcodng-Options.xlsx)
b) [Eden 2.5kb tiled amplicon protocol](https://www.pacb.com/wp-content/uploads/Customer-Collaboration-PacBio-Compatible-Eden-Protocol-for-SARS-CoV-2-Sequencing.pdf) and [PacBio barcoding options](https://www.pacb.com/wp-content/uploads/Eden-SARS-CoV-2-Primers-with-PacBio-Barcoding-Options.xlsx)
c) [CDC 500bp and 900bp tiled amplicon protocol](https://www.pacb.com/wp-content/uploads/Customer-Collaboration-PacBio-Compatible-CDC-Protocols-for-SARS-CoV-2-Sequencing.pdf) and [PacBio barcoding options](https://www.pacb.com/wp-content/uploads/CDC-SARS-CoV-2-Primers-with-PacBio-Barcoding-Options.xlsx)

Files are mirrored [here](../protocols/).

## _**4. IONTORRENT**_
A number of laboratories have reported success with the Ion AmpliSeq SARS-CoV-2 Research panel for the IonTorrent S5 platform. Amplicon strategies, such as ARTIC, should also work for the S5, and we'd welcome the addition of any working protocols and other resources to this section.
- [Ion AmpliSeq SARS-CoV-2 Research Panel for GeneStudio S5](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/MAN0019277_Ion_AmpliSeq_SARS-CoV-2_Research_Panel_GeneStudio_QR.pdf)

## _**5. TANDEM**_

### a) SARS-CoV-2 Parallel Sequencing by Illumina and ONT (UWMadison ZEST)
Staff and students from Thomas Friedrich and Dave O'Connor's laboratories at UWMadison have put together a tandem sequencing protocol and bioinformatic workflow that incorporates Illumina and ONT sequence. While this may be overkill for routine or high-throughput public health purposes, the necessary protocols, scripts and documentation are available here.
- [Protocol](https://openresearch.labkey.com/wiki/ZEST/Ncov/page.view?name=SARS0CoV-2%20Deep%20Sequencing)
- [Github resources](https://github.com/katarinabraun/SARS-CoV-2_sequencing)

## _**6. SANGER**_
### a) CDC NCIRD/DVD Sanger Tiling
An elegant approach from a more civilized age. The Viral Discovery laboratory team in CDC/NCIRD/DVD has used conventional Sanger sequencing to refine and complement betacoronavirus sequencing on next generation platforms.
- [Singleplex and Multiplex Protocols for Sanger/ONT/Illumina](./protocols/CDC-Comprehensive)
- [Preprint describing CDC protocols](https://www.biorxiv.org/content/10.1101/2020.04.22.055897v1)

---


# **Bioinformatic Tools, Scripts and Workflows**<a id="bioinformatics" />
### 0. General Resources.
A few different sites have resources that are foundational to most bioinformatics analysis.
- [NCBI](https://www.ncbi.nlm.nih.gov) and [EBI](https://www.ebi.ac.uk) are both INSDC public repositories, and contain all public access sequences, annotation, analysis and derived data. NCBI has added a specific landing page for SARS-CoV-2 research resources, available [here](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/).
- [GISAID](https://www.gisaid.org) was put into place to provide a framework for open sharing of influenza sequence data, while maintaining strict governance over the use of data and attribution to sequence submitters.
- The original SARS-CoV-2 reference sequence, "Wuhan-Hu-1", submitted by **Zhang Yongzhen and colleagues** is mirrored from NCBI in the `sequences` [folder](./sequences) here.  In that folder, you will also find a dated set of public genomes in that folder, downloaded from NCBI on a semi-regular basis. NB: These are consensus sequences manually pulled from [NCBI](https://www.ncbi.nlm.nih.gov) Genbank. More complete sets of sequences are available at the [GISAID](https://www.gisaid.org) public access repository, which requires account registration and adherence to a strict code of conduct.
- [Genexa](https://genexa.ch/sars2-bioinformatics-resources/) has put together a page with precomputed kmer sets, indicies, reference sequences, and a number of other useful resources for bioinformatic analysis of SARS-CoV-2 NGS data.
- [Illumina](https://www.illumina.com/informatics/specialized-bioinformatics-applications/coronavirus-software.html) has just released a SARS-CoV-2 Software Toolkit, including premade RT-PCR and sequencing workflows for ClarityLIMS, FPGA (DRAGEN)-powered bioinformatics, and direct submission to GISAID.

### 1. CDC NCIRD/DVD Bioinformatics SOPs
This section describes the basic bioinformatic workflow that the Viral Discovery laboratory in NCIRD, and other teams at CDC use for quality assessment, assembly and comparison of coronavirus sequences. IRMA, the Iterative Refinement Meta-Assembler developed by CDC's Influenza Division for routine influenza surveillance, has recently been updated to support both ebolavirus and coronavirus assembly tasks. While IRMA isn't used for all SARS-CoV-2 assemblies at CDC, it is a powerful tool for complex or problematic samples and datasets.
- [Preprint describing](https://www.biorxiv.org/content/10.1101/2020.04.22.055897v1)
- [IRMA: Iterative Refinement Meta-Assembler](https://wonder.cdc.gov/amd/flu/irma) is available here.
- CLCbio ([Jonathan Jacobs](mailto:jonathan.jacobs@qiagen.com)) has released workflows for the NCIRD/DVD protocols on [Illumina](./files/clc/CDC.Illumina.NEBNext.SARSCOV2.CLC.Workflow.v0.2.zip) and [Nanopore](./files/clc/CDC.Nanopore.Amplicon.SARSCOV2.Protocol.CLC.Worklow.v0.2.zip). Older versions available [here](./files/clc/).

### 2. CLCbio Genomics Workbench
QIAGEN has released example workflows and tutorials for analyzing Illumina and Oxford Nanopore SARS-CoV-2 sequence data using CLC Genomics Workbench v20.0.3. Note - these workflows are "Research Use Only" (RUO), and may need to modified to fit upstream protocols.  Free temporary licenses for CLC GWB and IPA are available, as well as a series of webinars and tutorials are available to familiarize users with the workflows. [Jonathan Jacobs](mailto:jonathan.jacobs@qiagen.com) and [Leif Schauser](mailto:leif.schauser@qiagen.com) are available for user support and specific questions. 

- [Temporary licenses for CLC Genomics Workbench and Ingenuity Pathway Analysis](https://go.qiagen.com/QDI-COVID19)
- [SARS-CoV-2 Analysis Tutorial with Nanopore Data](https://resources.qiagenbioinformatics.com/tutorials/Analysis_of_SARS-CoV-2.pdf) [(Figure)](https://user-images.githubusercontent.com/64888/77170622-ad637500-6a91-11ea-9b2a-f6f3cda4981e.png)
- [SARS-CoV-2 Workflow and Tutorial Data for Nanopore Data](http://resources.qiagenbioinformatics.com/testdata/SARS-CoV-2_MinION_example_data.zip) [(Figure)](https://user-images.githubusercontent.com/64888/77166238-4098ac80-6a8a-11ea-9b38-534fe956258b.png)
- [CLC Illumina Workflow v1](./files/clc/)
- [Additional CLC Genomics Workbench Workflows](https://github.com/jonathanjacobs/CLC-Genomics/tree/master/CLC-SARSCOV2-Workflows) developed by Jonathan Jacobs for QIAseq SARS-CoV-2 panel, Ion AmpliSeq SARS-COV-2 panel, and ARTICv3 Nanopore & llumina panels.

### 3. ARTIC Network Bioinformatics
The [ARTIC Network](https://artic.network) has released detailed instructions on how to setup and configure the conda environment needed to run their analysis pipelines. These are complete bioinformatic workflows, including runtime visualization, basecalling, mapping/assembly and reporting in a single, portable environment. The [artic-nCoV2019](https://github.com/artic-network/artic-ncov2019) repo includes source code and build instructions for a custom RAMPART configuration.  Additional instructions and documentation are available below.

- [Bioinformatics Environment Setup/Configuration](https://artic.network/ncov-2019/ncov2019-it-setup.html)

- [SARS-CoV-2 Bioinformatics SOP](https://artic.network/ncov-2019/ncov2019-bioinformatics-sop.html)

- [RAMPART Runtime Documentation](https://artic.network/ncov-2019/ncov2019-using-rampart.html)

- [Nextflow ARTIC nCoV Workflows - Connor Lab Cardiff/PHW](https://github.com/connor-lab/ncov2019-artic-nf)

### 4. BugSeq
[BugSeq](https://bugseq.com) has added support for automatic SARS-CoV-2 analysis (QC, consensus, variant calling and lineage detection) from nanopore sequencing data. This analysis ([example](https://bugseq.com/demo/sars-cov-2)) is triggered by detection of SARS-CoV-2 in data submissions and is automatically tailored for your experimental design (metagenomic and ARTIC v1,v2 and v3 amplicons supported). BugSeq is making analysis of SARS-CoV-2 samples available **free of charge** for any user publishing their data in a pre-print or the peer-reviewed literature. Additional information and documentation are available below:

- [Example report](https://bugseq.com/demo/sars-cov-2)

- [Additional documentation and bioinformatics details](https://docs.bugseq.com/tailored-experiment-analyses/covid-19-sars-cov-2)

### 5. One Codex
[One Codex](https://www.onecodex.com) has added support to its analysis platform for analyzing SARS-CoV-2 samples. This analysis ([example](https://app.onecodex.com/report/public/2b15c38901224e5f)) will be automatically run on any samples with SARS-CoV-2 reads. One Codex is making analysis of SARS-CoV-2 samples available **free of charge** for all users sharing their results and data publicly. Additional information and documentation are available below:

- [Blog post introducing COVID-19 analysis support](https://www.onecodex.com/blog/2020/03/16/covid-19-sequencing-analysis/)

- [An example report](https://app.onecodex.com/report/public/2b15c38901224e5f)

- [Additional documentation and bioinformatics details](https://docs.onecodex.com/en/articles/3793936-covid-19-sequencing-analysis)

### 6. Broad viral-ngs tools
The Broad Institute's viral genomics analysis tools can assist with assembly, metagenomics, QC, and NCBI submission prep, for Illumina-generated data on viral genomes. It is available in the following forms:

- The Terra cloud platform ([workspace](https://app.terra.bio/#workspaces/pathogen-genomic-surveillance/COVID-19) including example SARS-CoV-2 data from SRA, [blog post](https://support.terra.bio/hc/en-us/articles/360040613432), [getting started](https://support.terra.bio/hc/en-us/articles/360041068771)) 
- The DNAnexus cloud platform ([workflows](https://platform.dnanexus.com/projects/F8PQ6380xf5bK0Qk0YPjB17P/data/build/quay.io/broadinstitute/viral-pipelines/2.0.11.0-rc23))
- The Dockstore tool repository service - integrates with several cloud platforms, or download to run on-prem ([workflows](https://dockstore.org/organizations/BroadInstitute/collections/pgs))
- Github ([workflows](https://github.com/broadinstitute/viral-pipelines/tree/master/pipes/WDL))

The tools include:
- denovo and reference based assembly
- short read alignment and coverage plots
- krakenuniq metagenomic classification
- NCBI: SRA download, Genbank annotation download, Genbank submission prep
- multiple alignment of genomes w/MAFFT
- Illumina basecalling & demux, metrics, fastQC, ERCC spike-in counter

### 7. Genome Detective Virus tool

Genome Detective virus tool does QC, assembly and identification of SARS-CoV-2 from a wide range of sequencing protocols (metagenomic or targeted sequencing).

Raw sequence read files (FASTQ) can be uploaded directly in this [web-based tool](https://www.genomedetective.com/app/typingtool/virus/), and consensus sequences can be subsequently analyzed by the [the Coronavirus Typing Tool](https://www.genomedetective.com/app/typingtool/cov).

Example output:
- [Minion](https://www.genomedetective.com/app/typingtool/virus/job/example-SARSCoV2-nanopore)
- [Illumina](https://www.genomedetective.com/app/typingtool/virus/job/example-SARSCoV2-illumina)

### 8. CosmosID
[CosmosID](https://www.cosmosid.com) has recently posted a blog entry on their site, describing how to use their web-based analysis platform to analyze SARS-CoV-2 data.
- [Detection of SARS-CoV-2 Coronavirus using CosmosID](https://www.cosmosid.com/blog-in/2020/3/13/metagenomic-detection-of-sars-cov-2-coronavirus-using-cosmosid)

### 9. ARTIC on Illumina Bioinformatic Workflow
[@ErinYoung](https://github.com/erinyoung) and Kelly Oakeson at the Utah Department of Health have outlined their bioinformatics approach for SARS-CoV-2 sequences using ARTIC primers, sequenced on Illumina.
- [Bioinformatic protocol](./protocols/BFX-UT_ARTIC_Illumina)
- [iVar GitHub (Andersen Lab @ Scripps)](https://github.com/andersen-lab/ivar) and [documentation](https://andersen-lab.github.io/ivar/html/)

### 10. Galaxy
[Workflows](https://covid19.galaxyproject.org/) are available for the [Galaxy Platform](https://galaxyproject.org/)
- [Cheoinformatics: Virtual screening of the SARS-CoV-2 main protease](https://covid19.galaxyproject.org/cheminformatics/)
- [Genomics: Analysis of COVID-19 data using Galaxy, BioConda and public research infrastructure](https://covid19.galaxyproject.org/genomics)
- [ARTIC: Amplicon analysis using Artic workflows](https://covid19.galaxyproject.org/artic)

These repositories provide best practise workflows for genomic and chemoinformatic analyses for SARS-CoV-2 data. In addition to providing tools and workflows we provide free public computational infrastrcture for immediate use by anyone worldwide using a [consortium of Galaxy instances from US, EU, and Australia](https://galaxyproject.org/use/). 

The chemoinformatics workflows can be used to conduct fragment screening using molecular docking.  *This has [already been done for the SARS-COV-2 main protease (MPro)](https://covid19.galaxyproject.org/cheminformatics/histories/)* by the [Diamond Light Source's XChem team](https://www.diamond.ac.uk/Instruments/Mx/Fragment-Screening.html), [InformaticsMatters](http://informaticsmatters.com/) and the [European Galaxy Team](https://galaxyproject.eu/).  The genomics workflows use entirely open source software and open access platforms to perform (1) [data pre-processing](https://covid19.galaxyproject.org/genomics/1-PreProcessing), (2) [genome assembly](https://covid19.galaxyproject.org/genomics/2-Assembly), (3) [estimation of MRCA timing](https://covid19.galaxyproject.org/genomics/3-MRCA), (4) [analysis of intrahost variation](https://covid19.galaxyproject.org/genomics/4-Variation), (5) [analysis of substitutions within the *S* gene](https://covid19.galaxyproject.org/genomics/5-S-analysis), and (6) [analysis of recombination and selection](https://covid19.galaxyproject.org/genomics/6-RecombinationSelection). 

### 11. Nanopore Direct RNA Analysis Using MasterOfPores
The [Epitranscriptomics and RNA Dynamics Lab (Novoa)](https://www.crg.eu/en/programmes-groups/novoa-lab) and the Bioinformatics Core Facility ([BioCore](https://biocore.crg.eu/wiki/Main_Page) at the [CRG](https://www.crg.eu/) have released a set of tools and resources to support the analysis of nanopore direct RNA sequencing data.
- [Uniform analysis of coronavirus and SARS-COV-2 nanopore direct RNA sequencing datasets using MasterOfPores](https://biocorecrg.github.io/covid/)

### 12. BioNumerics SARS-CoV-2 Plugin
Applied Maths/bioMerieux have released a plugin for BioNumerics that facilitates the processing and analysis of SARS-CoV-2 genomic sequences, whether downloaded from a public data repository or generated locally. More info on the tool and a tutorial can be found [here](https://www.applied-maths.com/applications/sarscov2-genome-screening).

### 13. StaPH-B ToolKit Monroe Workflow
The [State Public Health Bioinformatics (StaPH-B)](http://www.staphb.org/) consortium has made their [Monroe workflow](https://staph-b.github.io/staphb_toolkit/workflow_docs/monroe/) accessible through the [StaPH-B ToolKit](https://staph-b.github.io/staphb_toolkit/). Monroe consists of three separate Nextflow pipelines for
- ARTIC + Illumina paired-end read assembly
- ARTIC + Oxford Nanopore Technlogies read assembly
- Cluster analysis from assembled SC2 genomes

### 14. fastv: identify SARS-CoV-2 from sequencing data in one minute
Fastv is a little-weight independent tool for ultra-fast identification of SARS-CoV-2 and other microbes from sequencing data. It detects SARS-CoV-2 sequences from FASTQ data, generates JSON reports and visualizes the result in HTML reports. It supports both short reads (Illumina, BGI, etc.) and long reads (ONT, PacBio, etc.). More information can be found [here](https://github.com/OpenGene/fastv).

### 15. EDGE COVID-19: A web platform for generating SARS-CoV-2 genomes

[EDGE COVID-19](https://edge-covid19.edgebioinformatics.org/) is a standardized web-based workflow for automated reference-based genome assembly of SARS-CoV-2 samples. The workflow accommodates Illumina or Oxford Nanopore Technologies data, performs read mapping and provides static and interactive figures/graphs to explore quality and any discovered SNPs, Variants, Gaps, and indels. Given raw FASTQ file(s) from amplicon-based methods (ARTIC, CDC) or shotgun sequencing (including from enrichment protocols), EDGE COVID-19 automates the production of a SARS-CoV-2 genome that is ready for submission to GISAID or GenBank. We have automated the process to submit high quality genomes to GISAID (with required metadata) and we are working on a similar process for SRA and GenBank.

Here is a [link](https://arxiv.org/abs/2006.08058) to the arXiv paper describing the workflow.

Our other COVID-19 related efforts can be found [here](https://covid19.edgebioinformatics.org/#/).

### 16. MiCall: Pipeline for processing NGS data to genotype human RNA viruses like SARS-CoV-2, HIV and hepatitis C
[MiCall] processes NGS read data from platforms like Illumina by either assembling them or mapping them to a set of
reference sequences. Then, it reports consensus sequences, variant mixtures, and quality control reports. For HIV and
hepatitis C, it also reports drug resistance interpretations of the variant mixtures.

MiCall is open-source software and comes packaged to be run under Docker or Singularity, to make installation easy. It's
also available to run on Illumina's [BaseSpace] web service, but the SARS-CoV-2 support is not yet available there.

[MiCall]: https://github.com/cfe-lab/MiCall
[BaseSpace]: https://basespace.illumina.com/apps/3247244?preferredVersion

---

# **Quality Management**<a id="quality" />
This section will describe best practices for laboratory and bioinformatic quality assurance, including preflight checks for sequence and metadata submission to public repositories.

### 1. FDA ARGOS SARS-CoV-2 Reference Sequence Data and Materials
The FDA, in collaboration with CDC, BEIR, UMaryland IGS, and others, have recently put together reference sequences, and materials for NGS sequencing and assay development.  These resources will be invaluable to many laboratories implementing NGS quality management programs. 
- [FDA ARGOS SARS-CoV-2 Resources](https://www.fda.gov/medical-devices/database-reference-grade-microbial-sequences-fda-argos/fda-argos-sars-cov-2-reference-grade-sequence-data)

### 2. Host Sequence Removal
SanitizeMe is a set of scripts and an X11 GUI for removing human host sequences from metagenomic data before SRA submission. Software and documentation are available [here](https://github.com/CDCgov/SanitizeMe).

- **to do**: I'd love for people to help describe their actual QC processes.

---


# **Submitting to Public Sequence Repositories**<a id="submitting" />

## WHO Code of Conduct: Sequence Sharing During Outbreaks
The draft *WHO code of conduct for open and timely sharing of pathogen genetic sequence data during outbreaks of infectious disease* lays out an important and sensible set of principles for sharing pathogen genetic sequence data during outbreaks of international importance.  The text is available [here](https://www.who.int/blueprint/what/norms-standards/GSDDraftCodeConduct_forpublicconsultation-v1.pdf?ua=1).

## Sequence naming conventions for public repositories
We are proposing simplified naming conventions for sequences submitted to GISAID and NCBI from US public health and clinical laboratories.

**COUNTRY** | /  **STATE-LAB-SAMPLE**  / | **YEAR**
--- | --- | ---
`USA` | ``/ CA-CDPH-S001 /`` | `2020`
`USA` | ``/ UT-UPHL-0601 /`` | `2020`
`USA` | ``/ AZ-TGEN-N1045 /`` | `2020`
`USA` | ``/ ID-UW-0316 /`` | `2020`

> The proposed convention is as follows: 1) country (USA). 2) The middle sample identification cell should include two-letter state (eg: CA), an abbreviated identifier for the submitting lab (eg: CDPH), as desired, and a unique sequence identifier (eg:01, S01, 454, ...), with all three terms separated by hyphens.




For states with only one submitting laboratory (which should be most), the identifier for the submitting laboratory may be omitted, resulting in a simple, state-level identifier such as `USA/UT-573/2020`.

*These recommendations are roughly compatible with existing submissions to GISAID and NCBI, but are completely open for debate.*
The current ICTV recommendations are [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7095448/pdf/41564_2020_Article_695.pdf), with the original biorXiv [here](https://www.biorxiv.org/content/10.1101/2020.02.07.937862v1.full.pdf).

## Recommended formatting and criteria for sample metadata

### [NCBI SARS-CoV-2 Genbank/SRA](https://www.ncbi.nlm.nih.gov/genbank/sars-cov-2-seqs/)
The National Center for Biotechnology has established a custom landing page for SARS-CoV-2 sequences and data, and is working to develop streamlined submission processes for Genbank and SRA.  For the time being, we suggest basing metadata and submission formatting on GISAID EpiCoV, which tends to be more comprehensive and structured. We will develop specific guidance for NCBI submissions.  In the meantime, here are some general resources to help with NCBI data submission and metadata management.

#### 1. NCBI Submission Portal
Individual sequences can be submitted to NCBI using the following web form. Create an NCBI user account, and select "SARS-CoV-2 (through BankIt)".
- [SARS-CoV-2 Submission Portal](https://submit.ncbi.nlm.nih.gov/sarscov2)

> NCBI has provided provisional guidance for SARS-CoV-2 sequence submissions to SRA and Genbank. [Detailed instructions are available here](./protocols/NCBI-Submission/NCBI-Sequence_submission_for_COVID-19.pdf). Any questions can be directed to NCBI staff [here](mailto:gb-admin@ncbi.nlm.nih.gov).

#### 2. NCBI Batch Submissions
NCBI has indicated that they plan to develop a specific rapid submission process for SARS-CoV-2 sequences. In the meantime, I believe you should be able to follow the FDA/CFSAN submission protocol below, which includes links to appropriate interfaces and templates (with obvious changes for pathogen and project information).

#### 3. FDA/CFSAN NCBI Submission and Data Curation Protocols
The FDA Center for Food Safety and Applied Nutrition ([CFSAN](https://www.fda.gov/about-fda/fda-organization/center-food-safety-and-applied-nutrition-cfsan)) has released a number of protocols as part of their [GenomeTrakr Network](https://www.fda.gov/food/whole-genome-sequencing-wgs-program/genometrakr-network) that may be useful for NCBI sequence submission and metadata curation. While they are written specifically for laboratories that are conducting routine sequencing of foodborne bacterial pathogens, these protocols provide an overview of sequence submission to the NCBI pathogen portal, metadata and preflight data checks.

- [NCBI submission protocol for microbial pathogen surveillance](https://www.protocols.io/view/ncbi-submission-protocol-for-microbial-pathogen-su-9aph2dn)
- [Populating the NCBI pathogen metadata template](https://www.protocols.io/view/populating-the-ncbi-pathogen-metadata-template-bck3iuyn)
- NCBI Data Curation - Pending release

### [GISAID EpiCoV](https://www.gisaid.org/)

The GISAID EpiCoV Public Access repository is based on existing submission processes and data structures for large-scale influenza surveillance (GISAID EpiFlu). As such, submitters to EpiCoV will discover that several of the required metadata submission fields may be problematic.  Nonetheless, a number of laboratories have been submitting sequences with the following:

METADATA FIELDS (GISAID) | GUIDANCE
--- | ---
`Virus name`| USA/FL-UF-103/2020 (see above)
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

A copy of the current bulk submission template is available [here](./templates/20200315_EpiCoV_BulkUpload_Template.xls).

---


# **Linking Sequence Accessions**<a id="linking" />
For data linkage, we are proposing the following template, as a simple, lightweight line list of tab-separated values. If this consensus recommendation for data linkage is acceptable, a preformatted .TSV will be made available. We recognize that not all samples sent for sequencing have a PUID associated.

SEQUENCE_NAME | GISAID_ID | GENBANK_ID | COLLECTION_DATE | PUID/COVID-ID
--- | --- | --- | --- | ---
USA/CA-CDPH-999/2020 | EPI_ISL_999999 | MT99999999 | 2020-04-01 | 99999

In this simple proposed schema, GISAID ID or GENBANK ID and COLLECTION DATE are required fields, and our hope is to maximize PUID completion.  All accession numbers, including PUID should be entered without any superfluous text or annotation.

---


# **Other Useful References and Resources**<a id="resources" />
## SARS-CoV-2 Genomics Consortia
- [UK COVID-19 Genomics Consortium (COG-UK)](https://cogconsortium.uk)

## Useful References
- [Open-source analytics tools for studying the COVID-19 coronavirus outbreak](https://www.medrxiv.org/content/10.1101/2020.02.25.20027433v2)
- [Global Initiative on Open-source Genomics for SARS-CoV-2](https://giogs.genomics.cn/)

## COVID-19 Host Genomics
- [COVID-19 Host Genetics Initiative](https://www.covid19hg.org)

## Slides and Presentations
- [StaPH-B SARS-CoV-2 Sequencing Seminar (Kevin Libuit, Virginia DCLS - 20200320)](./slides/SC2_Sequencing_and_Analysis_Libuit.pdf) [(Recording)](http://www.staphb.org/resources/staphb-videos/)

## Visualization and Phylogenetics
- [Nextstrain hCoV-19](https://www.nextstrain.org/ncov)
- [UCSC SARS-CoV-2 browser](https://genome.ucsc.edu/cgi-bin/hgTracks?db=wuhCor1)
- [CNCB 2019 Novel Coronavirus Resource (2019nCoVR)](https://bigd.big.ac.cn/ncov?lang=en)

## Diagnostic Resources
- [SARS-CoV-2 PCR Primer List (Grubaugh Lab, Yale SPH)](https://docs.google.com/spreadsheets/d/1kWLdaqWZZYt383EMSf-osyXGdeedrRF6YeBuOBin6Ok/edit#gid=0)
- [EDGE Bioinformatics: In-Silico Evaluation of Diagnostic Assays (LANL)](https://covid19.edgebioinformatics.org/#/assayValidation)

## Additional Protocol Resources
- [Coronavirus Methods Development Group @ protocols.io](https://www.protocols.io/groups/coronavirus-method-development-community)

---
# **Notices and Disclaimers**<a id="notices" />
## Public Domain

This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License

Unless otherwise specified, the repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

Any source code forked from other open source projects will inherit its license.

## Privacy

This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records

This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

---

Updated: 20200608 [@dmaccannell](https://www.twitter.com/dmaccannell)
