## Utah DoH ARTIC/Illumina Bioinformatic Workflow ([Erin Young](https://github.com/erinyoung)/Kelly Oakeson)
#### Uses ARTICv2 amplicons, sequenced on Illumina.

### Bioinformatic Workflow:
**STEP 1.** Map Illumina reads to the reference (MN908947.3; ARTIC default). This will also sort and remove unmapped reads:

```
bwa mem -t {threads} artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta {input.read1} {input.read2} | samtools sort | samtools view -F 4 -o {sample}.sorted.bam
```

**STEP 2.** Translate the ARTIC primer scheme into something readable by IVAR (perl wizard hat tip to [@torstenseemann](https://github.com/tseemann)):

```
perl -ne 'my @x=split m/\t/; print join("\t",@x[0..3], 60, $x[3]=~m/LEFT/?"+":"-"),"\n";' < artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.scheme.bed  > ARTIC-V1.bed
```

**STEP 3.** Trim the primers off of the bam sequences using ivar:

```
ivar trim -i {sample}.sorted.bam -b ARTIC-V1.bed -p {sample}.primertrim
```

**STEP 4.** Re-sort your bams:

```
samtools sort {sample}.primertrim -o {sample}.primertrim.sorted.bam
```

**STEP 5.** Get the consensus *fasta* that includes all the variants found, without replacing missing sequence with reference (missing sequence simply becomes "N"). The `samtools mpileup` options listed are those given in `ivar`'s [manual](https://andersen-lab.github.io/ivar/html/manualpage.html), and might not be the best options for our needs:

```
samtools mpileup -A -d 6000000 -B -Q 0 --reference artic-ncov2019/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta {sample}.primertrim.sorted.bam | ivar consensus -p {sample}.consensus -n N
```

**Possible QC options:**

Quast:
```
quast {sample}.consensus.fa -r artic-ncov2019/primer_schemes/nCoV-2019/V1/ncov-2019.reference.fasta --features GCF_009858895.2_ASM985889v3_genomic.gff --ref-bam {sample}.sorted.bam --output-dir quast/{sample}
```

samtools coverage (passing = > 95% Covered bases, Mean baseQ > 30, and Mean mapQ > 30):
```
samtools coverage {sample}.sorted.bam -o {sample}.samcov.txt
```
* side note: using the `-m` option prints out a text-based histogram that looks like:
```
MN908947.3 (29.9Kbp)
>  90.00% │▃███████████████████████████████ █████████████▇██▇│ Number of reads: 624800
>  80.00% │████████████████████████████████ █████████████████│ 
>  70.00% │████████████████████████████████ █████████████████│ Covered bases:   29.6Kbp
>  60.00% │████████████████████████████████ █████████████████│ Percent covered: 98.97%
>  50.00% │████████████████████████████████▆█████████████████│ Mean coverage:   4.49e+03x
>  40.00% │██████████████████████████████████████████████████│ Mean baseQ:      35.1
>  30.00% │██████████████████████████████████████████████████│ Mean mapQ:       60
>  20.00% │██████████████████████████████████████████████████│ 
>  10.00% │██████████████████████████████████████████████████│ Histo bin width: 598bp
>   0.00% │██████████████████████████████████████████████████│ Histo max bin:   100%
          1        6.0K     12.0K     17.9K     23.9K      29.9K  
```

### iVar Analysis Package (Andersen Lab @ Scripps)
- [iVar GitHub Repository](https://github.com/andersen-lab/ivar)
- [iVar Documentation](https://andersen-lab.github.io/ivar/html/)
