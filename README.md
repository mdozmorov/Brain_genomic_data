# Brain-related -omics data

## `Allen_Brain_Atlas` - Allen Brain Atlas, Human brain

RNA-Sequencing datasets from two brains

http://human.brain-map.org/static/download

These datasets contain gene expression values (raw and TPM counts) for a selected set of anatomic structures matched across the two brains, as well as sample and gene metadata necessary for analysis: H0351.2001, H0351.2002.

- `AllenReferenceAtlas_v1_2008_102011.pdf` - abbreviations of anatomical structures, http://help.brain-map.org/download/attachments/2818169/AllenReferenceAtlas_v1_2008_102011.pdf?version=1

Brain regions:

- PL	parietal lobe
- FL	frontal lobe
- OL	occipital lobe
- Ins	insula
- TL	temporal lobe
- PHG	parahippocampal gyrus
- Str	striatum
- CgG	cingulate gyrus
- GP	globus pallidus
- CbCx	cerebellar cortex

- `Contents.txt` - describes the files in the zipped folders.

- `rnaseq_donor9861.zip` - 121 samples. `RNAseqTPM.csv` - first column - gene name, `SampleAnnot.csv`

```
CbCx  CgG   FL   GP  Ins   OL  PHG   PL  Str   TL 
   4    4   39    4    4    9    4   28    9   16 
```

- `rnaseq_donor10021.zip` - 121 samples

```
CbCx  CgG   FL   GP  Ins   OL  PHG   PL  Str   TL 
   5    4   37    4    4    8    4   31    8   16 
```


## `BrainSpan` - BrainSpan transcriptomics and epigenomics of the developing human brain

Brain regions:

- DFC	Dorsolateral prefrontal cortex
- VFC	Ventrolateral prefrontal cortex
- MFC	Anterior (rostral) cingulate (medial prefrontal) cortex
- OFC	Orbital frontal cortex
- M1C	Primary motor cortex (area M1, area 4)
- S1C	Primary somatosensory cortex (area S1, areas 3,1,2)
- IPC	Posteroinferior (ventral) parietal cortex
- A1C	Primary auditory cortex (core)
- STC	Posterior (caudal) superior temporal cortex (area TAc)
- ITC	Inferolateral temporal cortex (area TEv, area 20)
- V1C	Primary visual cortex (striate cortex, area V1/17)
- HIP	Hippocampus (hippocampal formation)
- AMY	Amygdaloid complex
- STR	Striatum
- MD	Mediodorsal nucleus of thalamus
- CBC	Cerebellar cortex
- M1C-S1C	Primary motor-sensory cortex (samples)

### `Methylation`

http://download.alleninstitute.org/brainspan/Methylation/ 

- `Methylation_ReadMe.doc` - brief description, 450K arrays
- `1109_methylation_beta_values.zip` - 207M, `wget http://download.alleninstitute.org/brainspan/Methylation/1109_methylation_beta_values.zip`. 485592 rows x 87 columns. Sample IDs like "7796806001_R01C01.AVG_BETA"
- `1110_methylation_beta_values.zip` - 219M, `wget http://download.alleninstitute.org/brainspan/Methylation/1110_methylation_beta_values.zip`

Unzipped files are tab-delimited txt files, approx. twice in size, with commented header followed by the matrix of cgIDs and beta values. 

### `microRNA`

http://download.alleninstitute.org/brainspan/MicroRNA/

- `MicroRNA.xls` - Raw RNA-seq counts, 1861 miRNAs in 216 samples, http://download.alleninstitute.org/brainspan/MicroRNA/MicroRNA.xls 

### `RNA-seq`

This data set contains RNA-Seq RPKM (reads per kilobase per million; see the whitepaper at www.brainspan.org) values averaged to genes.

- `genes_matrix_csv.zip` - 62M, RNA-Seq Gencode v10 summarized to genes, http://www.brainspan.org/api/v2/well_known_file_download/267666525. compressed folder with files: 

- columns_metadata.csv -- the samples are listed in the same order as the columns in expression_matrix. 524 samples, age, gender, structure

```
    A1C     AMY      CB     CBC     CGE     DFC     DTH     HIP     IPC     ITC     LGE     M1C M1C-S1C      MD 
     31      33       3      29       2      35       5      32      33      34       2      26       5      24 
    MFC     MGE     Ocx     OFC     PCx     S1C     STC     STR     TCx     URL     V1C     VFC 
     32       2       2      31       2      26      36      28       1       2      33      35 
```

- rows_metadata.csv -- the genes are listed in the same order as the rows in expression_matrix.csv. 52376 genes, ensembl_gene_id, gene_symbol, entrez_id

- expression_matrix.csv -- the rows are genes and the columns samples; the first column is the row number. 52376 genes in 524 samples


## `CommonMind`

https://www.synapse.org/#!Synapse:syn2759792/wiki/

Controlled access data. >1000 postmortem brain samples from donors with Schizophrenia, Bipolar disease and individuals with no neuropsychiatric disorders - originating from tissue collections at four brain banks. Data consists of DNA and RNA sequencing, genotyping and epigenetics. 	Dorsolateral Prefrontal Cortex.

- [Data Terms of Use](https://www.synapse.org/#!Synapse:syn2759792/wiki/197282)


## `Roadmap`

Histone and DNAse hypersensitive sites, BED and WIG files, http://genboree.org/EdaccData/Release-9/sample-experiment/

Brain regions:

- Cortex derived primary cultured neurospheres
- Ganglion Eminence derived primary cultured neurospheres
- Brain Angular Gyrus
- Brain Anterior Caudate
- Brain Cingulate Gyrus
- Brain Germinal Matrix
- Brain Hippocampus Middle
- Brain Inferior Temporal Lobe
- Brain_Dorsolateral_Prefrontal_Cortex
- Brain Substantia Nigra
- Fetal Brain Male
- Fetal Brain Female


## `PsychENCODE`

https://www.synapse.org//#!Synapse:syn4921369/wiki/235539

Controlled access data. Regulatory genomic elements (promoters, enhancers, silencers and insulators), catalog epigenetic modifications and quantify coding and non-coding RNA and protein expression in tissue and cell-type-specific samples from healthy (neurotypical) control, disease-affected post-mortem human brains, and complementary investigations of induced pluripotent cells (iPSCs), cultured neuronal cells derived from olfactory neuroepithelium (CNON cells), and model organisms. The project currently focuses on three major psychiatric disorders: Autism Spectrum Disorder (ASD), Bipolar Disorder (BD) and Schizophrenia (SCZ), prioritizing brain regions and cell types that previous research has suggested contribute to these disorders.

- RNA-, miRNA-, lncRNA-, ChIP-, SNP-, ATAC-, WGBS-, HiC- seq datasets, https://www.synapse.org/#!Synapse:syn4921369/wiki/390659

## `HumanFC`

https://www.synapse.org/#!Synapse:syn5321694

ATAC-seq data - BAM, BED, BigWig, Metadata, Controlled access data. ATAC-seq has been generated on the DLPFC from 84 control, 82 Schizophrenia and 4 Bipolar Disorder adult brain samples.

## `EpiMap`

https://www.synapse.org/#!Synapse:syn4566010

ChIP-seq on H3K27ac and H3K4me3, Controlled access data. Mostly normal brains, PFC and ACC from the same individual.


## `FANTOM5` 

http://slidebase.binf.ku.dk/human_enhancers/presets

BED files of enhancers in general "brain" tissue and in specific cells, e.g. neuron, dendritic cell, astrocyte.


## `Misc`

A recent inventory of 2,104 human- gained enhancers active during cerebral corticogenesis [@Reilly:2015aa]. https://www.ncbi.nlm.nih.gov/pubmed/25745175, ChIP-seq of H3K27ac, H3K4me2, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63649
