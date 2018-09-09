# Brain-related omics data

Mostly human, but some datasets contain omics data from model organisms.

* [Consortia](#consortia)
  * [Allen Brain Atlas](#allen-brain-atlas)
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
* [Hi-C data](#hi-c-data)
* [Single-cell RNA-seq](#single-cell-rna-seq)
* [Misc](#misc)

# Consortia

## Allen Brain Atlas

Allen Brain Atlas, Human brain. RNA-Sequencing datasets from two brains. http://human.brain-map.org/static/download

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


## BrainSpan

BrainSpan transcriptomics and epigenomics of the developing human brain

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

The Brain eQTL Almanac. Gene- or SNP-centric exploration of brain region-specific eQTLs

http://www.braineac.org/

## BrainGVEX

Genetic variants affect brain gene expression and risks of psychiatric disorders. eQTLs, pQTLs, csQTLs (chromatin state) from genotypes, RNA-seq, ATAC-seq, protein arrays

https://www.synapse.org/#!Synapse:syn4590909, Controlled access data 

# Hi-C data

- Two human brain regions, dorsolateral prefrontal cortex, hippocampus. No replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87112
- Two human brain regions: the cortical and subcortical plate (CP), consisting primarily of post-mitotic neurons and the germinal zone (GZ), containing primarily mitotically active neural progenitors. Three replicates per condition. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77565.
- Six human brain tumours: five glioblastomas ( GB176, GB180, GB182, GB183 and GB238) and one anaplastic astrocytoma (AA86). No replicates. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE81879
- Mouse neural progenitors (NPCs), and cortical neurons (CNs), purified NPC and CN populations from neocortex (ncx_NPC, ncx_CN). Replicates (4-6). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96107.

# Single-cell RNA-seq

- Darmanis, S., Sloan, S.A., Zhang, Y., Enge, M., Caneda, C., Shuer, L.M., Hayden Gephart, M.G., Barres, B.A., and Quake, S.R. (2015). A survey of human brain transcriptome diversity at the single cell level. Proc. Natl. Acad. Sci. USA 112, 7285–7290. - Single cell brain transcriptomics, human. Fluidigm C1 platform. Healthy cortex cells (466 cells) containing: Astrocytes, oligodendrocytes, oligodendrocyte precursor cells (OPCs), neurons, microglia, and vascular cells. Single cells clustered into 10 clusters, their top 20 gene signatures are in Supplementary Table S3. Raw data athttps://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67835
    - `data/TableS3.txt` - top 20 cell type-specific genes
    - `data/TableS3_matrix.txt` - genes vs. cell types with 0/1 indicator variables.

- Zeisel, A., Munoz-Manchado, A.B., Codeluppi, S., Lonnerberg, P., La Manno, G., Jureus, A., Marques, S., Munguba, H., He, L., Betsholtz, C., et al. (2015). Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA- seq. Science 347, 1138–1142. - 3,005 single cells from the hippocampus and cerebral cortex of mice. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361, http://linnarssonlab.org/cortex/, and more on this site.
    - `data/Zeisel_2015_TableS1.xlsx` - Table S1 - gene signatures for Ependymal, Oligodendrocyte, Microglia, CA1 Pyramidal, Interneuron, Endothelial, S1 Pyramidal, Astrocyte, Mural cells. [Source](http://science.sciencemag.org/highwire/filestream/628248/field_highwire_adjunct_files/1/aaa1934_TableS1.xlsx)

- The 1 million neuron data set from the 10X Genomics website (hhttps://support.10xgenomics.com/single-cell-gene-expression/datasets). R packages for its analyses: `TENxGenomics`, https://github.com/mtmorgan/TENxGenomics, `TENxBrainAnalysis`, https://github.com/Bioconductor/TENxBrainAnalysis

- Luo, Chongyuan, Christopher L. Keown, Laurie Kurihara, Jingtian Zhou, Yupeng He, Junhao Li, Rosa Castanon, et al. “Single-Cell Methylomes Identify Neuronal Subtypes and Regulatory Elements in Mammalian Cortex.” Science (New York, N.Y.) 357, no. 6351 (11 2017): 600–604. https://doi.org/10.1126/science.aan3351. - single-cell methylation of human and mouse neuronal cells. Marker genes with cell type-specific methylation profiles - Table S3, http://science.sciencemag.org/content/suppl/2017/08/09/357.6351.600.DC1

- Nowakowski, Tomasz J., Aparna Bhaduri, Alex A. Pollen, Beatriz Alvarado, Mohammed A. Mostajo-Radji, Elizabeth Di Lullo, Maximilian Haeussler, et al. “Spatiotemporal Gene Expression Trajectories Reveal Developmental Hierarchies of the Human Cortex.” Science 358, no. 6368 (December 8, 2017): 1318–23. https://doi.org/10.1126/science.aap8809. - single-cell RNA-seq of neuronal cell types. Dimensionality reduction, clustering, WGCNA, defining cell type-specific signatures, comparison with other signatures (Zeng, Miller). 
    - `data/Nowakowski_2017_Tables_S1-S11.xlsx` - Table S5 has brain region-specific gene signatures. [Source](http://science.sciencemag.org/highwire/filestream/703290/field_highwire_adjunct_files/1/aap8809_Nowakowski_SM-Tables-S1-S11.xlsx) 

- Major Depressive Disorder Working Group of the Psychiatric Genomics Consortium et al., “Genetic Identification of Brain Cell Types Underlying Schizophrenia,” Nature Genetics 50, no. 6 (June 2018): 825–33, https://doi.org/10.1038/s41588-018-0129-5. - Cell-type specificity of schizophrenia SNPs judged by enrichment in expressed genes. scRNA-seq custom data collection. Difference between schizophrenia and neurological disorders.
    - `data/Brain_cell_type_gene_expression.xlsx` - Supplementary Table 4 - Specificity values for Karolinska scRNA-seq superset. Specificity represents the proportion of the total expression of a gene found in one cell type as compared to that in all cell types (i.e., the mean expression in one cell type divided by the mean expression in all cell types). Gene X cell type matrix. Level 1 (core cell types) and level 2 (extended collection of cell types) data. [Source](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0129-5/MediaObjects/41588_2018_129_MOESM3_ESM.xlsx)

# Misc

- `data/Brain_genes_TableS1.xlsx` - mouse brain cell type-specific genes. From Budakian et al. - 2014 - Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq. [Source](http://www.sciencemag.org/content/suppl/2015/02/18/science.aaa1934.DC1/aaa1934_TableS1.xlsx)

- A recent inventory of 2,104 human-gained enhancers active during cerebral corticogenesis [@Reilly:2015aa]. https://www.ncbi.nlm.nih.gov/pubmed/25745175, ChIP-seq of H3K27ac, H3K4me2, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63649

- Strand, Andrew D., Aaron K. Aragaki, Zachary C. Baquet, Angela Hodges, Philip Cunningham, Peter Holmans, Kevin R. Jones, Lesley Jones, Charles Kooperberg, and James M. Olson. “Conservation of Regional Gene Expression in Mouse and Human Brain.” PLoS Genetics 3, no. 4 (April 20, 2007): e59. https://doi.org/10.1371/journal.pgen.0030059. - Normal human and mouse brain gene expression comparison, _motor cortex_, _caudate nucleus_, _cerebellum_ for human and analogous anterior cortex, striatum, and cerebellum samples. Microarrays, 12 human samples in 3 conditions in one set (GSE3790) and 9 human samples in 2 conditions in another set. 6 mouse samples in 3 conditions. Supplement: Normalized Affy probe ID data for individual subjects (human, mouse). Differentially expressed genes between each pair of brain regions http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0030059#s5. Human normal brain samples are from a large Huntington's disease study, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3790

- Mo, Alisa, Eran A. Mukamel, Fred P. Davis, Chongyuan Luo, Gilbert L. Henry, Serge Picard, Mark A. Urich, et al. “Epigenomic Signatures of Neuronal Diversity in the Mammalian Brain.” Neuron 86, no. 6 (June 17, 2015): 1369–84. https://doi.org/10.1016/j.neuron.2015.05.018. - mouse brain, three neuronal populations. RNA-seq, ATAC-seq, histone ChIP-seq. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63137

- The Brainstorm Consortium, Verneri Anttila, Brendan Bulik-Sullivan, Hilary K. Finucane, Raymond K. Walters, Jose Bras, Laramie Duncan, et al. “Analysis of Shared Heritability in Common Disorders of the Brain.” Science 360, no. 6395 (June 22, 2018): eaap8757. https://doi.org/10.1126/science.aap8757. - Heritability correlation among psychiatric and neurological disorders, contrasted with phenotypes. Neurologic and psychiatric disorders are distinct. Details about each pair of correlations. Table S7 - all pairwise correlations, Table S13 - data sources. http://science.sciencemag.org/content/suppl/2018/06/20/360.6395.eaap8757.DC1
    - `data/Brainstorm_Table_S7.xlsx` - Table S7. Disorder-disorder (A), disorder-phenotype (B) and phenotype-phenotype (C) correlation results. [Source](http://science.sciencemag.org/highwire/filestream/711735/field_highwire_adjunct_files/4/aap8757_Table_S7.xlsx)
    - `data/Brainstorm_Table_S13.xlsx` - Table S13. Data sources, responsible consortia, and data availability. [Source](http://science.sciencemag.org/highwire/filestream/711735/field_highwire_adjunct_files/3/aap8757_Table_S13.xlsx)



Issues and/or Pull requests to add other data are welcome!
