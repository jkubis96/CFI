## CFI: Cell Functionality and Interaction Analysis Tool

![Python version](https://img.shields.io/badge/python-%E2%89%A53.12-blue?logo=python&logoColor=white.png) ![License](https://img.shields.io/badge/license-GPLv3-blue) ![Docs](https://img.shields.io/badge/docs-available-blueviolet)




<p align="right">
<img  src="https://github.com/jkubis96/Logos/blob/main/logos/jbs_current.png?raw=true" alt="drawing" width="200" />
</p>


### Author: Jakub Kubiś 

<div align="left">
 Institute of Bioorganic Chemistry<br />
 Polish Academy of Sciences<br />
</div>


## Description




</br>


CFI (Cell Functionality and Interaction Analysis Tool) is an analytical platform designed for the investigation of cellular functions, intra-cellular gene interactions, and intercellular communication, including ligand–receptor interactions and adhesion junctions.

The tool integrates the advanced enrichment and interaction analysis capabilities of [GEDSpy](https://github.com/jkubis96/GEDSpy) with the single-cell data processing functionalities provided by [JDtI](https://github.com/jkubis96/JDtI).

CFI extends these capabilities by enabling the identification of direct cell–cell interactions, dominant biological processes, and gene–gene interaction networks within individual cells, along with comprehensive visualization of these relationships.


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/log.png" alt="drawing" width="500" />
</p>


CFI is designed to streamline the interpretation of biological data, enabling researchers to perform in-depth analyses of cellular functions and interactions for drug target discovery.


</br>

Included data bases:

* [Gene Ontology (GO-TERM)](http://geneontology.org/)
* [KEGG (Kyoto Encyclopedia of Genes and Genomes)](https://www.genome.jp/kegg/)
* [Reactome](https://reactome.org/)
* [HPA (Human Protein Atlas)](https://www.proteinatlas.org/)
* [NCBI](https://www.ncbi.nlm.nih.gov/)
* [STRING](https://string-db.org/)
* [IntAct](https://www.ebi.ac.uk/intact/home)
* [CellTalk](https://tcm.zju.edu.cn/celltalkdb/)
* [CellPhone](https://www.cellphonedb.org/)


*If you use CFI, please remember to cite both CFI and the original sources of the data you utilized in your work.*

In the case of interactions network analyses, it is recommended to use the  [JVectorGraph](https://github.com/jkubis96/JVectorGraph) library to easily adjust and customise graph and network visualisations from the Python side.


<br />

## 📚 Table of Contents

- [Installation](#installation)
- [Documentation](#doc)
- [Example usage](#example)
  - [1. Cell functionality](#bf)
    - [1.1. Create project](#bf1)
    - [1.2. Calculate cell marker genes](#bf2)
    - [1.3. Marker gene enrichment analysis](#bf3)
    - [1.4. Cell inside gene interactions](#bf4)
    - [1.5. Cell-cell interactions](#bf5)
    - [1.6. Saving & loading project](#bf6)
  - [2. Comparison of cell interaction sets](#br)
    - [2.1. Create projects](#br1)
    - [2.2. Calculate interactions](#br2)
    - [2.3. Comparison analysis](#br3)


<br />


# Installation <a id="installation"></a>

#### In command line write:

```
pip install cfi-toolkit
```




## Documentation <a id="doc"></a>


Documentation for classes and functions is available here 👉 [Documentation 📄](https://jkubis96.github.io/CFI/cfi.html)


<br />


## Example usage <a id="example"></a>

### 1. Cell functionality <a id="bf"></a>

##### 1.1. Create project <a id="bf1"></a>

```
# ------------------------------------------------------------
# Import required standard and project-specific libraries
# ------------------------------------------------------------

import os
from jdti import COMPsc          # JDtI module for handling single-cell projects
from cfi_toolkit import CellFunCon       # Cell functional connectivity / enrichment analysis


# ------------------------------------------------------------
# Load single-cell sequencing data using the JDtI framework
# ------------------------------------------------------------

# Define the project directory containing input data
# and specify the sample identifiers to be loaded

jseq_object = COMPsc.project_dir(
    os.path.join(os.getcwd(), "data"),   # path to data directory
    ["s1"]                               # list of project/sample IDs
)

# Load sparse expression matrices from the project structure
# normalized_data=True ensures that pre-normalized counts are used

jseq_object.load_sparse_from_projects(normalized_data=True)


# ------------------------------------------------------------
# Initialize CellFunCon analysis object
# ------------------------------------------------------------

# Create a CellFunCon instance using the loaded single-cell dataset.
# This object enables downstream functional enrichment,
# interaction analysis, and pathway inference.

instance = CellFunCon(jseq_object)
```

<br />


##### 1.2. Calculate cell marker genes <a id="bf2"></a>

```
# ------------------------------------------------------------
# Required step before functional enrichment analysis
# ------------------------------------------------------------
# Identify cell-type–specific marker genes based on:
# - minimum expression threshold (min_exp)
# - minimum fraction of expressing cells (min_pct)
# - parallel computation across multiple processes (n_proc)
#
# The resulting marker set is used as input for downstream
# functional enrichment and pathway analysis.

instance.calculate_cells_markers(
    min_exp=0,
    min_pct=0.05,
    n_proc=10
)
```

<br />


##### 1.3. Marker gene enrichment analysis <a id="bf3"></a>

```
# ------------------------------------------------------------
# Perform functional enrichment analysis for each cell population
# ------------------------------------------------------------
# This step evaluates overrepresentation of biological functions,
# pathways, and gene sets using previously identified marker genes.
#
# Parameters:
# - p_value : significance threshold for enrichment results
# - log_fc  : minimum log fold-change required for marker inclusion
# - top_max : maximum number of top-ranked markers used per cell type
#
# The output provides functionally annotated cell states
# for downstream biological interpretation.

instance.enrich_cells_fucntionality(
    p_value=0.05,
    log_fc=0.25,
    top_max=500
)
```

<br />


* GO-TERM

```
data = instance.get_enrichment_data( 
                        data_type = 'GO-TERM', 
                        p_value = 0.05, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)
```

| parent     | parent_genes                                   | parent_pval_FISH | parent_pval_BIN | parent_n | parent_pct | child     | child_genes                              | child_pval_FISH | child_pval_BIN | child_n | child_pct | parent_name                         | child_name                                              | cell             | source  |
|------------|-------------------------------------------------|------------------|-----------------|----------|------------|-----------|-------------------------------------------|------------------|----------------|---------|-----------|-------------------------------------|----------------------------------------------------------|------------------|---------|
| GO:0060090 | ['LETMD1', 'TNRC6A']                            | 0.01397          | 0.01402         | 2        | 0.04348    | GO:0030674 | ['TNRC6A', 'SRRT', 'CRADD', 'SLMAP']      | 3.22e-05         | 3.30e-05       | 4       | 0.08696   | MF : molecular adaptor activity     | MF : protein-macromolecule adaptor activity              | STRIATUM_1 # s1 | GO-TERM |
| GO:0032991 | ['CLU','SRRT','HSP90AA1','SNX4','CHEK1','PPP3CB'] | 9.38e-05         | 9.49e-05        | 6        | 0.13043    | GO:0005955 | ['PPP3CB']                                | 0.00459          | 0.00458        | 1       | 0.02174   | CC : protein-containing complex     | CC : calcineurin complex                                | STRIATUM_1 # s1 | GO-TERM |
| GO:0032991 | ['CLU','SRRT','HSP90AA1','SNX4','CHEK1','PPP3CB'] | 9.38e-05         | 9.49e-05        | 6        | 0.13043    | GO:0031428 | ['FBL']                                   | 0.00535          | 0.00535        | 1       | 0.02174   | CC : protein-containing complex     | CC : box C/D methylation guide snoRNP complex           | STRIATUM_1 # s1 | GO-TERM |
| GO:0032991 | ['CLU','SRRT','HSP90AA1','SNX4','CHEK1','PPP3CB'] | 9.38e-05         | 9.49e-05        | 6        | 0.13043    | GO:0005664 | ['ORC6']                                  | 0.00611          | 0.00611        | 1       | 0.02174   | CC : protein-containing complex     | CC : nuclear origin of replication recognition complex  | STRIATUM_1 # s1 | GO-TERM |
| GO:0032991 | ['CLU','SRRT','HSP90AA1','SNX4','CHEK1','PPP3CB'] | 9.38e-05         | 9.49e-05        | 6        | 0.13043    | GO:0008287 | ['PPP3CB']                                | 0.00611          | 0.00611        | 1       | 0.02174   | CC : protein-containing complex     | CC : protein serine/threonine phosphatase complex       | STRIATUM_1 # s1 | GO-TERM |

<br />


 Visualization

```
from cfi_toolkit import encrichment_cell_heatmap

fig = encrichment_cell_heatmap(data = data,
                             fig_size = (3,3), 
                             sets = None,
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/heatmap.svg" alt="drawing" width="450" />
</p>


<br />

* KEGG

```
data = instance.get_enrichment_data( 
                        data_type = 'KEGG', 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)
```
| 2nd                     | 2nd_genes                                                                 | 2nd_pval | 2nd_n | 2nd_pct | 3rd                     | 3rd_genes                                   | 3rd_pval | 3rd_n | 3rd_pct | cell             | source |
|-------------------------|---------------------------------------------------------------------------|----------|-------|---------|--------------------------|-----------------------------------------------|----------|-------|---------|------------------|--------|
| Neurodegenerative disease | ['ATXN10','VDAC2','PSMA3','APC','SETX','CREBBP','PSMD8','PSMD7','NDUFA4'] | 1.10e-05 | 9     | 0.093   | Spinocerebellar ataxia   | ['ATXN10','VDAC2','PSMA3','PSMD8','PSMD7']    | 7.51e-06 | 5     | 0.052   | STRIATUM_2 # s1 | KEGG   |
| Neurodegenerative disease | ['ATXN10','VDAC2','PSMA3','APC','SETX','CREBBP','PSMD8','PSMD7','NDUFA4'] | 1.10e-05 | 9     | 0.093   | Huntington disease       | ['VDAC2','PSMA3','CREBBP','PSMD8','PSMD7','NDUFA4'] | 2.08e-05 | 6 | 0.062 | STRIATUM_2 # s1 | KEGG |
| Neurodegenerative disease | ['ATXN10','VDAC2','PSMA3','APC','SETX','CREBBP','PSMD8','PSMD7','NDUFA4'] | 1.10e-05 | 9     | 0.093   | Alzheimer disease        | ['VDAC2','PSMA3','APC','PSMD8','PSMD7','NDUFA4'] | 8.20e-05 | 6 | 0.062 | STRIATUM_2 # s1 | KEGG |
| Neurodegenerative disease | ['ATXN10','VDAC2','PSMA3','APC','SETX','CREBBP','PSMD8','PSMD7','NDUFA4'] | 1.10e-05 | 9     | 0.093   | Prion disease            | ['VDAC2','PSMA3','PSMD8','PSMD7','NDUFA4']    | 1.23e-04 | 5     | 0.052   | STRIATUM_2 # s1 | KEGG   |
| Neurodegenerative disease | ['ATXN10','VDAC2','PSMA3','APC','SETX','CREBBP','PSMD8','PSMD7','NDUFA4'] | 1.10e-05 | 9     | 0.093   | Parkinson disease        | ['VDAC2','PSMA3','PSMD8','PSMD7','NDUFA4']    | 1.77e-04 | 5     | 0.052   | STRIATUM_2 # s1 | KEGG   |


<br />

Visualization

```
from cfi_toolkit import encrichment_cell_heatmap

fig = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = None,
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/heatmap2.svg" alt="drawing" width="450" />
</p>


> [!NOTE]
> In this case, **no significant KEGG terms** were detected for the second cell type.  
> To indicate the absence of specific terms, use  
> `instance1.get_included_cells()` to display both cells at the top.

<br />

Visualization

```
from cfi_toolkit import encrichment_cell_heatmap

fig2 = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)
```

<br />

 Visualization



<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/heatmap3.svg" alt="drawing" width="450" />
</p>


<br />


* REACTOME

```
data = instance.get_enrichment_data( 
                        data_type = 'REACTOME', 
                        p_value = 0.05, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)
```

| pathway                                      | genes        | p-value | n | pct  | top level pathway        | top genes            | cell             | source   |
|----------------------------------------------|-------------|--------|---|------|---------------------------|----------------------|------------------|----------|
| Resistance of ERBB2 KD mutants to tesevatinib | [HSP90AA1]  | 0.00306 | 1 | 0.0217 | Disease                   | [RNGTT,HSP90AA1,HMGB1] | STRIATUM_1 # s1 | REACTOME |
| Resistance of ERBB2 KD mutants to osimertinib | [HSP90AA1]  | 0.00306 | 1 | 0.0217 | Disease                   | [RNGTT,HSP90AA1,HMGB1] | STRIATUM_1 # s1 | REACTOME |
| Resistance of ERBB2 KD mutants to AEE788      | [HSP90AA1]  | 0.00306 | 1 | 0.0217 | Disease                   | [RNGTT,HSP90AA1,HMGB1] | STRIATUM_1 # s1 | REACTOME |
| Apoptosis induced DNA fragmentation           | [HMGB1]     | 0.00991 | 1 | 0.0217 | Programmed Cell Death     | [HSP90AA1,HMGB1]       | STRIATUM_1 # s1 | REACTOME |
| Regulation of CDH11 mRNA translation by miRNAs| [TNRC6A]    | 0.00839 | 1 | 0.0217 | Cell-Cell communication   | [TNRC6A]               | STRIATUM_1 # s1 | REACTOME |

<br />

Visualization

```
from cfi_toolkit import encrichment_cell_heatmap

fig = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/heatmap4.svg" alt="drawing" width="450" />
</p>


<br />

* Specificity (HPA)

```
data = instance.get_enrichment_data( 
                        data_type = 'specificity', 
                        p_value = 1, 
                        test = 'FISH', 
                        adj = 'BH', 
                        parent_inc = False, 
                        top_n = 50)
```
| specificity            | genes                | p-value | n | pct  | cell             | source      |
|------------------------|----------------------|--------|---|------|------------------|-------------|
| Cardiomyocytes         | [SLMAP]              | 0.116  | 1 | 0.0217 | STRIATUM_1 # s1 | specificity |
| Proximal tubular cells | [ALDH6A1]            | 0.098  | 1 | 0.0217 | STRIATUM_1 # s1 | specificity |
| granulocytes           | [CLU, ALDH6A1, TIMP2]| 0.132  | 3 | 0.0652 | STRIATUM_1 # s1 | specificity |
| liver                  | [ALDH6A1]            | 0.288  | 1 | 0.0217 | STRIATUM_1 # s1 | specificity |
| kidney                 | [ALDH6A1]            | 0.146  | 1 | 0.0217 | STRIATUM_1 # s1 | specificity |
| Schwann cells          | [GRIK3]              | 0.030  | 1 | 0.0217 | STRIATUM_1 # s1 | specificity |

<br />

Visualization

```
from cfi_toolkit import encrichment_cell_heatmap

fig = encrichment_cell_heatmap(data = enr,
                             fig_size = (3,3), 
                             sets = instance1.get_included_cells(),
                             top_n = 3,
                             test = 'FISH', 
                             adj = 'BH', 
                             parent_inc = False,
                             font_size = 16,
                             clustering = 'ward',
                             scale = True)
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/heatmap5.svg" alt="drawing" width="450" />
</p>

<br />

##### 1.4. Cell inside gene interactions <a id="bf4"></a>

View available cell names

```
# ------------------------------------------------------------
# List all cell populations included in the dataset
# ------------------------------------------------------------
# Useful for verifying dataset structure before selecting
# specific cell populations for further investigation.

instance.get_included_cells()
```

<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/out_cell.bmp" alt="drawing" width="450" />
</p>



Retrieve interaction data for the selected cell

```
cell_int = instance.get_gene_interactions('STRIATUM_1 # s1')
```

| A    | B       | interaction_type     | connection_type     | source      |
|------|---------|----------------------|---------------------|-------------|
| CLU  | CLU     | physical association | protein -> protein  | Alzheimers  |
| CLU  | HSP90AA1| physical association | protein -> protein  | Alzheimers  |
| FBL  | KRR1    |                      | gene -> gene        | STRING      |
| FBL  | KRR1    |                      | protein -> protein  | STRING      |
| KRR1 | FBL     |                      | gene -> gene        | STRING      |
| KRR1 | FBL     |                      | protein -> protein  | STRING      |


<br />

 Visualization

```
from cfi_toolkit import gene_interaction_network

fig5 = gene_interaction_network(idata = cell_int, min_con = 2)

from JVG import JVG

nt = JVG.NxEditor(fig5)
nt.edit()
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/gin.bmp" alt="drawing" width="300" />
</p>

<br />


##### 1.5. Cell-cell interactions <a id="bf5"></a>


Calculate cells interactions

```
# ------------------------------------------------------------
# Infer functional and molecular connections between cell populations
# ------------------------------------------------------------
# This step identifies potential interactions based on
# inferred molecular communication signals
#
# The resulting network describes relationships between cells
# and can be used for downstream visualization, pathway analysis,
# or biological interpretation of tissue organization.

instance.calculate_cell_connections()
```

Get data

```
cell_con = instance.get_cell_connections()
```

| interaction              | directionality     | classification                 | modulatory_effect | interactor1 | interactor2 | cell1       | cell2       |
|--------------------------|--------------------|---------------------------------|-------------------|-------------|-------------|-------------|-------------|
| CDH2 -> CDH2             | Adhesion-Adhesion  | Adhesion by Cadherin           |                   | CDH2        | CDH2        | STRIATUM_1  | STRIATUM_2  |
| CDH6 -> CDH6             | Adhesion-Adhesion  | Adhesion by Cadherin           |                   | CDH6        | CDH6        | STRIATUM_1  | STRIATUM_2  |
| CDH7 -> CDH7             | Adhesion-Adhesion  | Adhesion by Cadherin           |                   | CDH7        | CDH7        | STRIATUM_1  | STRIATUM_2  |
| COL11A1 -> ITGA1+ITGB1   | Adhesion-Adhesion  | Adhesion by Collagen/Integrin  |                   | COL11A1     | ITGA1       | STRIATUM_1  | STRIATUM_2  |
| COL11A1 -> ITGA1+ITGB1   | Adhesion-Adhesion  | Adhesion by Collagen/Integrin  |                   | COL11A1     | ITGB1       | STRIATUM_1  | STRIATUM_2  |

<br />

Visualization

```
from cfi_toolkit import draw_cell_conections

fig = draw_cell_conections(cell_con)

from JVG import JVG

nt = JVG.NxEditor(fig)
nt.edit()
```

<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/cell_con.bmp" alt="drawing" width="300" />
</p>

<br />



##### 1.6. Saving & loading project <a id="bf6"></a>


Saving current project

```
instance.save_project('project')
```


Loading previously saved project

```
from cfi_toolkit import CellFunCon   

instance = CellFunCon.load_project('project.psc')
```


<br />


### 2. Comparison of cell interaction sets <a id="br"></a>

##### 2.1. Create projects <a id="br1"></a>

```
# ------------------------------------------------------------
# Import required libraries
# ------------------------------------------------------------
import os
from jdti import COMPsc          # JDtI module for single-cell project handling
from cfi_toolkit import CellFunCon       # Functional analysis and cell interaction inference


# ------------------------------------------------------------
# Load single-cell datasets for two experimental conditions
# ------------------------------------------------------------
# Each COMPsc object represents an independent project/sample
# loaded from the same data directory but identified by
# different sample IDs (e.g., control vs. case).

jseq_object1 = COMPsc.project_dir(
    os.path.join(os.getcwd(), "data"),
    ["s1"]   # first dataset / condition
)
jseq_object1.load_sparse_from_projects(normalized_data=True)


jseq_object2 = COMPsc.project_dir(
    os.path.join(os.getcwd(), "data"),
    ["s2"]   # second dataset / condition
)
jseq_object2.load_sparse_from_projects(normalized_data=True)


# ------------------------------------------------------------
# Initialize CellFunCon analysis objects for comparative study
# ------------------------------------------------------------
# Separate instances enable:
# - independent marker detection
# - condition-specific functional enrichment
# - downstream comparison of cell functionality and interactions

instance1 = CellFunCon(jseq_object1)
instance2 = CellFunCon(jseq_object2)
```
<br />


##### 2.2. Calculate interactions <a id="br2"></a>

```
# ------------------------------------------------------------
# Calculate cell–cell connection networks for each dataset
# ------------------------------------------------------------
# This step reconstructs potential functional and molecular
# interactions between cell populations independently for:
# - dataset 1 (e.g., control)
# - dataset 2 (e.g., disease / treated condition)

instance1.calculate_cell_connections()
instance2.calculate_cell_connections()
```

<br />


##### 2.3. Comparison analysis <a id="br3"></a>

```
# ------------------------------------------------------------
# Compare cell–cell interaction networks between conditions
# ------------------------------------------------------------
# This step performs a comparative analysis of inferred
# cell–cell connections across multiple biological conditions.
#
# Each entry in `instances_dict` represents a fully processed
# CellFunCon object with reconstructed interaction networks.
# Here, we contrast:
# - healthy condition
# - disease condition

from cfi_toolkit import compare_connections

instances_dict = {
    "healthy": instance2,
    "disease": instance1
}

comparison = compare_connections(instances_dict=instances_dict, 
                                 cells_compartment = None, 
                                 connection_type  = ['Adhesion-Adhesion',
                                                     'Gap-Gap',
                                                     'Ligand-Ligand',
                                                     'Ligand-Receptor',
                                                     'Receptor-Receptor',
                                                     'Undefined'])
```

| feature | p-value  | pct_valid | pct_ctrl | FC    | log(FC) | norm_diff | group   |
|---------|----------|-----------|----------|-------|---------|-----------|---------|
| HSP90B1 | 4.75e-14 | 0.164     | 0.935    | 0.160 | -2.65   | -6.58     | healthy |
| CANX    | 5.90e-14 | 0.091     | 0.848    | 0.102 | -3.30   | -5.53     | healthy |
| ARPC5   | 1.34e-11 | 0.055     | 0.739    | 0.080 | -3.64   | -4.93     | healthy |
| CALR    | 5.47e-09 | 0.182     | 0.761    | 0.226 | -2.14   | -4.41     | healthy |
| APLP2   | 1.16e-08 | 0.073     | 0.609    | 0.115 | -3.12   | -3.73     | healthy |
| ITGB1   | 1.44e-08 | 0.055     | 0.587    | 0.098 | -3.36   | -3.61     | healthy |

<br />

Visualization

```
from cfi_toolkit import volcano_plot_conections

fig = volcano_plot_conections(
    deg_data = comparison,
    p_adj = True,
    top = 25,
    top_rank = "p_value",
    p_val = 0.05,
    lfc = 0.25,
    rescale_adj = True,
    image_width = 12,
    image_high = 12,
)
```


<p align="center">
<img  src="https://github.com/jkubis96/CFI/raw/refs/heads/lib_merging/fig/volcano.svg" alt="drawing" width="450" />
</p>


<br />

An example analysis pipeline is available here → [Example file](example.py)

<br />


#### Have fun JBS©







