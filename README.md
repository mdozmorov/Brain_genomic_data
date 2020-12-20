# Brain-related omics data

Mostly human, but some datasets contain omics data from model organisms.

* [Consortia](#consortia)
  * [Allen Brain Atlas](#allen-brain-atlas)
    * [Human Brain Atlas](#human-brain-atlas)
    * [Aging, dementia and TBI](#aging--dementia-and-tbi)
    * [IVY Glioblastoma Atlas](#ivy-glioblastoma-atlas)
  * [BrainSpan](#brainspan)
    * [Methylation](#methylation)
    * [microRNA](#microrna)
    * [RNA-seq](#rna-seq)
  * [CommonMind](#commonmind)
  * [Roadmap](#roadmap)
  * [PsychENCODE](#psychencode)
  * [HumanFC](#humanfc)
  * [EpiMap](#epimap)
  * [FANTOM5](#fantom5)
  * [Braineac](#braineac)
  * [BrainGVEX](#braingvex)
* [ATAC-seq data](#atac-seq-data)
* [Hi-C data](#hi-c-data)
* [scRNA-seq data](#scrna-seq-data)
* [Misc](#misc)

# Consortia

## Allen Brain Atlas

- All Allen Brain Atlases and data, human and mouse. http://portal.brain-map.org
    - **Cell types database** (electrophysiological, morphological, and transcriptomic data measured from individual cells)
    - **Brain observatory** (in vivo recordings)
    - **Mouse brain connectivity atlas** (image data detailing cell types in different brain areas)
    - **Reference atlases** (anatomy maps of human and mouse brain)
    - **Mouse brain atlas**, **Developing mouse brain atlas** (IHC gene expression, spatial correlation)
    - **Adult and developing non-human primate (NHP) atlas** (gene expression, neuroanatomical data)
    - **Mouse spinal cord atlas** (gene expression and imaging data)

### Human brain atlas

- Hawrylycz, Michael J., Ed S. Lein, Angela L. Guillozet-Bongaarts, Elaine H. Shen, Lydia Ng, Jeremy A. Miller, Louie N. van de Lagemaat, et al. “[An Anatomically Comprehensive Atlas of the Adult Human Brain Transcriptome](https://doi.org/10.1038/nature11405).” Nature 489, no. 7416 (September 2012) - [Allen Human Brain Atlas](https://portal.brain-map.org/) - gene expression of approx. 900 neuroanatomical slices of brain from two individuals. In Situ Hybridization, Microarrays, MRI. Gene expression correlate with spatial localization. [Data download](http://human.brain-map.org/static/)download
    - Microarray datasets. All normalized microarray expression values as well as probe and sample metadata.
    - RNA-Sequencing datasets from two brains. Gene expression values (raw and TPM counts) for a selected set of anatomic structures matched across the two brains, as well as sample and gene metadata.
    - [AllenReferenceAtlas_v1_2008_102011.pdf](http://help.brain-map.org/download/attachments/2818169/AllenReferenceAtlas_v1_2008_102011.pdf?version=1) - abbreviations of anatomical structures. Brain regions: **PL** - parietal lobe; **FL** - frontal lobe; **OL** occipital lobe; **Ins** - insula; **TL** - temporal lobe; **PHG** - parahippocampal gyrus; **Str** striatum; **CgG** cingulate gyrus; **GP** - globus pallidus; **CbCx** cerebellar cortex.
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

### Aging, dementia and TBI

The Aging, Dementia and Traumatic Brain Injury Study is a detailed neuropathologic, molecular and transcriptomic characterization of brains of control and TBI exposure cases from a unique aged population-based cohort from the Adult Changes in Thought (ACT) study, http://aging.brain-map.org/

De-identified clinical information, protein quantification, normalized and unnormalized RNA-seq FPKM values. http://aging.brain-map.org/download/index

### IVY Glioblastoma Atlas

The Ivy Glioblastoma Atlas Project is a foundational resource for exploring the anatomic and genetic basis of glioblastoma at the cellular and molecular levels. http://glioblastoma.alleninstitute.org/

Patient and tumor information, normalized and unnormalized RNA-seq FPKM values. http://glioblastoma.alleninstitute.org/static/download.html


## BrainSpan

BrainSpan transcriptomics and epigenomics of the developing human brain. http://www.brainspan.org/

Brain regions: **DFC** - Dorsolateral prefrontal cortex; **VFC** - Ventrolateral prefrontal cortex; **MFC** - Anterior (rostral) cingulate (medial prefrontal) cortex; **OFC** - Orbital frontal cortex; **M1C** - Primary motor cortex (area M1, area 4); **S1C** - Primary somatosensory cortex (area S1, areas 3,1,2); **IPC** - Posteroinferior (ventral) parietal cortex; **A1C** - Primary auditory cortex (core); **STC** - Posterior (caudal) superior temporal cortex (area TAc); **ITC** - Inferolateral temporal cortex (area TEv, area 20); **V1C** - Primary visual cortex (striate cortex, area V1/17); **HIP** - Hippocampus (hippocampal formation); **AMY** - Amygdaloid complex; **STR** - Striatum; **MD** - Mediodorsal nucleus of thalamus; **CBC** - Cerebellar cortex; **M1C-S1C** - Primary motor-sensory cortex (samples)

### Methylation

http://download.alleninstitute.org/brainspan/Methylation/ 

- `Methylation_ReadMe.doc` - brief description, 450K arrays
- `1109_methylation_beta_values.zip` - 207M, `wget http://download.alleninstitute.org/brainspan/Methylation/1109_methylation_beta_values.zip`. 485592 rows x 87 columns. Sample IDs like "7796806001_R01C01.AVG_BETA"
- `1110_methylation_beta_values.zip` - 219M, `wget http://download.alleninstitute.org/brainspan/Methylation/1110_methylation_beta_values.zip`

Unzipped files are tab-delimited txt files, approx. twice in size, with commented header followed by the matrix of cgIDs and beta values. 

### microRNA

http://download.alleninstitute.org/brainspan/MicroRNA/

- `MicroRNA.xls` - Raw RNA-seq counts, 1861 miRNAs in 216 samples, http://download.alleninstitute.org/brainspan/MicroRNA/MicroRNA.xls 

### RNA-seq

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


## CommonMind

https://www.synapse.org/#!Synapse:syn2759792/wiki/

Controlled access data. >1000 postmortem brain samples from donors with Schizophrenia, Bipolar disease and individuals with no neuropsychiatric disorders - originating from tissue collections at four brain banks. Data consists of DNA and RNA sequencing, genotyping and epigenetics. 	Dorsolateral Prefrontal Cortex.

- [Data Terms of Use](https://www.synapse.org/#!Synapse:syn2759792/wiki/197282)


## Roadmap

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


## PsychENCODE

https://www.synapse.org//#!Synapse:syn4921369/wiki/235539

Controlled access data. Regulatory genomic elements (promoters, enhancers, silencers and insulators), catalog epigenetic modifications and quantify coding and non-coding RNA and protein expression in tissue and cell-type-specific samples from healthy (neurotypical) control, disease-affected post-mortem human brains, and complementary investigations of induced pluripotent cells (iPSCs), cultured neuronal cells derived from olfactory neuroepithelium (CNON cells), and model organisms. The project currently focuses on three major psychiatric disorders: Autism Spectrum Disorder (ASD), Bipolar Disorder (BD) and Schizophrenia (SCZ), prioritizing brain regions and cell types that previous research has suggested contribute to these disorders.

- RNA-, miRNA-, lncRNA-, ChIP-, SNP-, ATAC-, WGBS-, HiC- seq datasets, https://www.synapse.org/#!Synapse:syn4921369/wiki/390659

- PsychENCODE resource. Genomics of the human brain: 1866 samples, gene expression (bulk and ~32,000 single-cell, decomposed with NMF into cell type profiles), enhancers, Hi-C, eQTLs and eGenes, GWAS variants linked to genes, regulatory networks, integrated data using deep learning. Samples from healthy and psychiatric disorder patients. Four data levels: raw, processed, derived, integrative. Download: gene expression matrix (bulk and single-cell), differentially expressed and spliced (cross-disease) genes, enhancer (H3K27ac) peaks, QTL maps. http://resource.psychencode.org
    - Wang, Daifeng, Shuang Liu, Jonathan Warrell, Hyejung Won, Xu Shi, Fabio C. P. Navarro, Declan Clarke, et al. “Comprehensive Functional Genomic Resource and Integrative Model for the Human Brain.” Science 362, no. 6420 (December 14, 2018): eaat8464. https://doi.org/10.1126/science.aat8464.

- Human brain development, integrative analysis by PsychEncode and BrainSpan. Genotyping, bulk and scRNA-seq, small RNA-seq, histone-seq, methylation. Detailed analysis of pre- and postnatal periods differences, figures. Processed and raw data at http://development.psychencode.org/ and http://brainspan.org/
    - Li, Mingfeng, Gabriel Santpere, Yuka Imamura Kawasawa, Oleg V. Evgrafov, Forrest O. Gulden, Sirisha Pochareddy, Susan M. Sunkin, et al. “Integrative Functional Genomic Analysis of Human Brain Development and Neuropsychiatric Risks.” Science 362, no. 6420 (December 14, 2018): eaat7615. https://doi.org/10.1126/science.aat7615.


## HumanFC

https://www.synapse.org/#!Synapse:syn5321694

ATAC-seq data - BAM, BED, BigWig, Metadata, Controlled access data. ATAC-seq has been generated on the DLPFC from 84 control, 82 Schizophrenia and 4 Bipolar Disorder adult brain samples.

## EpiMap

https://www.synapse.org/#!Synapse:syn4566010

ChIP-seq on H3K27ac and H3K4me3, Controlled access data. Mostly normal brains, PFC and ACC from the same individual.

## FANTOM5

http://slidebase.binf.ku.dk/human_enhancers/presets

BED files of enhancers in general "brain" tissue and in specific cells, e.g. neuron, dendritic cell, astrocyte.

## Braineac

The Brain eQTL Almanac. Gene- or SNP-centric exploration of brain region-specific eQTLs. 130 individuals with multiple brain region samples. The aim of Braineac is to release to the scientific community a valid instrument to investigate the genes and SNPs associated with neurological disorders. http://www.braineac.org/

## BrainGVEX

Genetic variants affect brain gene expression and risks of psychiatric disorders. eQTLs, pQTLs, csQTLs (chromatin state) from genotypes, RNA-seq, ATAC-seq, protein arrays

https://www.synapse.org/#!Synapse:syn4590909, Controlled access data 

# ATAC-seq data

- [MouseBrain](http://catlas.org/mousebrain/#!/home) - scATAC-seq from 45 regions of mouse brain. Data explorer, text-based download. http://catlas.org/mousebrain/#!/home
    - Li, Yang Eric, Sebastian Preissl, Xiaomeng Hou, Ziyang Zhang, Kai Zhang, Rongxin Fang, Yunjiang Qiu, et al. “An Atlas of Gene Regulatory Elements in Adult Mouse Cerebrum.” Preprint. Neuroscience, May 11, 2020. https://doi.org/10.1101/2020.05.10.087585.
 

# Hi-C data

- Two human brain regions, dorsolateral prefrontal cortex, hippocampus. No replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112
- Two human brain regions: the cortical and subcortical plate (CP), consisting primarily of post-mitotic neurons and the germinal zone (GZ), containing primarily mitotically active neural progenitors. Three replicates per condition. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77565.
- Six human brain tumours: five glioblastomas ( GB176, GB180, GB182, GB183 and GB238) and one anaplastic astrocytoma (AA86). No replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81879
- Mouse neural progenitors (NPCs), and cortical neurons (CNs), purified NPC and CN populations from neocortex (ncx_NPC, ncx_CN). Replicates (4-6). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107.

# scRNA-seq data

See [scRNA-seq notes, Brain section](https://github.com/mdozmorov/scRNA-seq_notes#brain)

- `data/Brain_scRNA-seq_TableS3.txt` - brain cell type-specific genes, from Darmanis S. et.al, and Stephen R. Quake. “A Survey of Human Brain Transcriptome Diversity at the Single Cell Level.” PNAS 2015 https://www.ncbi.nlm.nih.gov/pubmed/26060301. Ten signatures, 20-genes each. `data/Brain_scRNA-seq_TableS3_matrix.txt` - signature matrix reformatted into genes vs. cell types, with each cell having 1/0 to indicate a gene is a part or not of a cell type-specific signature.

- `data/Brain_genes_TableS1.xlsx` - mouse brain cell type-specific genes. From Budakian et al. - 2014 - Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. [Source](http://www.sciencemag.org/content/suppl/2015/02/18/science.aaa1934.DC1/aaa1934_TableS1.xlsx)



# Misc

- [NeMo](https://nemoarchive.org/) - the Neuroscience Multi-omics Data Archive. Data organized by assay (chromatin, methylation, transcriptome), grant, lab, organism (human, mouse). [Direct download](http://data.nemoarchive.org/)

- Epigenomics of brain metastasis - Illumina 450K methylation of 96 brain metastasis specimens from patients with lung, breast, and melanoma cancers. BrainMETH - a classifier to distinguish molecular subtypes. minfi data processing. [Raw and normalized data at GSE108576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108576)
    - Salomon et al., “[Brain Metastasis DNA Methylomes, a Novel Resource for the Identification of Biological and Clinical Features](https://doi.org/10.1038/sdata.2018.245)”

- meQTLs using whole genome bisulfite sequencing in schizophrenia, DLPFC and hippocampus regions. 38% CpGs are meQTLs. Expermiental confounders have little effect. Accounting for LD, clustering into regions. Results from two regions overlap. Comparison with controls, DMRs. Gene ontology enrichment of meQTL-associated genes in brain-related processes. CpH DNAm level is less associated with genomic variations.Age associations. All results are in the [supplementary tables](https://www.biorxiv.org/content/10.1101/2020.09.24.311878v1.supplementary-material)
    - Perzel Mandell, Kira A, Nicholas J Eagles, Richard Wilton, Amanda J Price, Stephen A Semick, Leonardo Collado-Torres, Ran Tao, et al. “[Widespread Methylation Quantitative Trait Loci and Their Role in Schizophrenia Risk](https://doi.org/10.1101/2020.09.24.311878).” Preprint. Genomics, September 24, 2020.

- `Neurodata.io` - Neurodata.io is a repository of image data. R packages and tools for statistical data analysis. https://neurodata.io/tools/
    - Burns, Randal, Eric Perlman, Alex Baden, William Gray Roncal, Ben Falk, Vikram Chandrashekhar, Forrest Collman, et al. “A Community-Developed Open-Source Computational Ecosystem for Big Neuro Data.” ArXiv:1804.02835 [q-Bio], April 9, 2018. http://arxiv.org/abs/1804.02835.

- A recent inventory of 2,104 human-gained enhancers active during cerebral corticogenesis [@Reilly:2015aa]. https://www.ncbi.nlm.nih.gov/pubmed/25745175, ChIP-seq of H3K27ac, H3K4me2, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63649

- Strand, Andrew D., Aaron K. Aragaki, Zachary C. Baquet, Angela Hodges, Philip Cunningham, Peter Holmans, Kevin R. Jones, Lesley Jones, Charles Kooperberg, and James M. Olson. “Conservation of Regional Gene Expression in Mouse and Human Brain.” PLoS Genetics 3, no. 4 (April 20, 2007): e59. https://doi.org/10.1371/journal.pgen.0030059. - Normal human and mouse brain gene expression comparison, _motor cortex_, _caudate nucleus_, _cerebellum_ for human and analogous anterior cortex, striatum, and cerebellum samples. Microarrays, 12 human samples in 3 conditions in one set (GSE3790) and 9 human samples in 2 conditions in another set. 6 mouse samples in 3 conditions. Supplement: Normalized Affy probe ID data for individual subjects (human, mouse). Differentially expressed genes between each pair of brain regions http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030059#s5. Human normal brain samples are from a large Huntington's disease study, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3790

- Mo, Alisa, Eran A. Mukamel, Fred P. Davis, Chongyuan Luo, Gilbert L. Henry, Serge Picard, Mark A. Urich, et al. “Epigenomic Signatures of Neuronal Diversity in the Mammalian Brain.” Neuron 86, no. 6 (June 17, 2015): 1369–84. https://doi.org/10.1016/j.neuron.2015.05.018. - mouse brain, three neuronal populations. RNA-seq, ATAC-seq, histone ChIP-seq. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63137

- The Brainstorm Consortium, Verneri Anttila, Brendan Bulik-Sullivan, Hilary K. Finucane, Raymond K. Walters, Jose Bras, Laramie Duncan, et al. “Analysis of Shared Heritability in Common Disorders of the Brain.” Science 360, no. 6395 (June 22, 2018): eaap8757. https://doi.org/10.1126/science.aap8757. - Heritability correlation among psychiatric and neurological disorders, contrasted with phenotypes. Neurologic and psychiatric disorders are distinct. Details about each pair of correlations. Table S7 - all pairwise correlations, Table S13 - data sources. http://science.sciencemag.org/content/suppl/2018/06/20/360.6395.eaap8757.DC1
    - `data/Brainstorm_Table_S7.xlsx` - Table S7. Disorder-disorder (A), disorder-phenotype (B) and phenotype-phenotype (C) correlation results. [Source](http://science.sciencemag.org/highwire/filestream/711735/field_highwire_adjunct_files/4/aap8757_Table_S7.xlsx)
    - `data/Brainstorm_Table_S13.xlsx` - Table S13. Data sources, responsible consortia, and data availability. [Source](http://science.sciencemag.org/highwire/filestream/711735/field_highwire_adjunct_files/3/aap8757_Table_S13.xlsx)



Issues and/or Pull requests to add other data are welcome!
