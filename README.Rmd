# DoReMiTraExplorer

**DoReMiTra-explorer** is an interactive Shiny application for exploring and visualizing radiation transcriptomic datasets from the DoReMiTra collection. It enables users to intuitively investigate transcriptional responses to radiation through a suite of customizable plots, including PCA, dose–response gene plots, boxplots, and heatmaps of highly variable genes.

## Installation

``` r
# Install the DoReMiTra package

BiocManager::install("DoReMiTra")

# Install the DoReMiTra-explorer Shiny App

Clone the repository:

git clone https://github.com/AhmedSAHassan/DoReMiTra-shiny.git
```

## Key Features

### **Principal Component Analysis** (PCA)

- Explore global transcriptomic structure

- Use custom number of top variable genes

- Color samples by metadata such as Dose, Time_point, or Sex, etc,.

### **Boxplots** to compare expression distributions across conditions

### Dose–response **gene plots** to examine radiation-induced trends

- allow selection of multiple genes of interest

- Color samples by metadata

- Links to the respective GeneCards page is provided for the selected genes

### **Heatmaps** of Highly Variable Genes

User-adjustable number of genes

Optional hierarchical clustering


## Example usage

``` r
library(DoReMiTra)

# Load an example SummarizedExperiment from DoReMiTra

se <- get_DoReMiTra_data("SE_Salah_2025_ExVivo")


# Launch the app

DoReMiTra_explorer(se)
```


## Links

DoReMiTra R Package (GitHub):
https://github.com/AhmedSAHassan/DoReMiTra

