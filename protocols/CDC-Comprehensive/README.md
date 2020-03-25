# CDC Comprehensive SARS-CoV-2 Sequencing Protocols

Thie is a comprenesive document including the protocols used by the DVD Pathogen Discovery lab at CDC for sequencing full genomes of SARS-CoV-2. Details are given for both a singleplex, tiling, nested RT-PCR that can be used with Sanger, ONT, or Illumina sequencing, as well as a multiplex, tiling RT-PCR that can be used with ONT or Illumina sequencing.  

Additionally, there are support files for primer trimming, pirmer ordering, and consensus generation.

Further support files will come along soon.
## Files
- Multplex_PCR/ - Requred primers for Multiplex PCR protocol
    + SC2_200324.bedpe - BED PE file with primer locations, used fro trimming
    + SC2_200324.fasta - Primers for 6-pool Multiplex PCR.
    + SC2_200324.tsv
- Singleplex_PCR/ - Required for Singplex PCR and Sanger sequencing
    +  Additional_Seq_primers.fasta - Additional primers required for Sanger sequencing
    +  Additional_Seq_primers.tsv 
    +  Singleplex_Primers.fasta - Primers for nested, singleplex PCR
    +  Singleplex_Primers.tsv
- scripts/ - Files to help with data processing
- 
## Table of Contents
1. Disclaimers
1. Singleplex nested RT-PCR
    1. Protocol Notes
    1. Required Reagents
    1. Procedure
        1. First round of RT-PCR
        1. Second round of semi-nested or nested PCR
1. Sanger sequencing
    1. Required Reagents
    1. Procedure
1. Multiplex PCR
    1. Protocol Notes
    1. Required reagents
    1. Procedure
        1. Generate primer pools
        1. First-strand synthesis
        1. Multiplex PCR
1. Nanopore Sequencing
    1. Protocol Notes
    1. Procedure
        1. Barcode amplicons
        1. Prepare Nanopore Ligation-based Library
        1. Load MinION and sequence
        1. Generate consensus sequences from MinION data
        1. Quality control and analysis suggestions
1. Illumina Library Preparation and Sequencing
    1. Protocol Notes
    1. Required Reagents
    1. Procedure for Library Preparation
        1.	Fragmentation and End Repair
        1.	Adapter Ligation
        1.	PCR enrichment of Adapter-Ligated DNA
        1.	Sizing and quantitation
1. MiSeq sequencing
    1. Protocol Notes
    1. Required Reagents
    1. Procedure
        1. Dilute and Pool Libraries
        1. Denature Libraries
        1. Load and Run MiSeq
        1. Generation of consensus sequences from MiSeq data
1. Appendix A – Singleplex PCR Primers
1. Appendix B – Sequencing Primers
1. Appendix C – Plate Setup for Nested PCR and Sanger Sequencing
1. Appendix D – Multiplex PCR Primers
1. Appendix E – AMPure XP bead clean-up
1. Appendix F – Quantitation using Qubit
    1. Required reagents
    1. Procedure
1. Appendix G – CENTRI-SEP 96 Protocol

