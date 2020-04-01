# CDC Comprehensive SARS-CoV-2 Sequencing Protocols

This is a comprehensive set documents, including the protocols used by the DVD Pathogen Discovery lab at CDC for sequencing full genomes of SARS-CoV-2. Details are given for both a single-plex tiling nested RT-PCR that can be used with Sanger, ONT or Illumina sequencing platforms, as well as a multiplex tiling RT-PCR that can be used with either ONT or Illumina sequencing.  

Additionally, there are support files for primer trimming, pirmer ordering, and consensus sequence generation.

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

- [CDC_SARS-CoV-2_Sequencing_200325-2.pdf](./CDC_SARS-CoV-2_Sequencing_200325-2.pdf) - complete documented protocols (March 25, 2020).
