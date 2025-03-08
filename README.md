# CNS_CTan

### Overview

This project is a comprehensive attempt to annotate a kidney single-cell RNA sequencing (scRNA-seq) dataset using CellTypist, a machine learning-based tool for automatic cell type annotation. The dataset, obtained from the HuBMAP consortium, was originally intended to be processed and annotated using Azimuth (R-based); however, due to various technical and time constraints, the analysis was instead completed using Python.

This README.md provides an in-depth walkthrough of the steps taken, challenges faced, and the reasoning behind critical decisions made throughout the project.

### Learning About CellTypist and AnnData

Before starting the implementation, I invested time in understanding CellTypist and AnnData by exploring official documentation, research articles, and online resources. Some of the key learnings include:

AnnData (Annotated Data): A widely used data structure in Python for storing and managing single-cell omics data. It efficiently handles gene expression matrices along with associated metadata.

CellTypist: A deep learning-based cell type annotation tool that leverages pre-trained models to predict cell identities from single-cell transcriptomic data.

Gene Symbol Mapping: Understanding the importance of converting Ensembl Gene IDs to HGNC symbols to ensure compatibility with CellTypist models.

By leveraging online tutorials, GitHub repositories, and documentation, I familiarized myself with these tools to apply them effectively to the dataset.

### Detailed Step-by-Step Walkthrough of the Project

This section provides an extremely detailed breakdown of the steps I followed to complete the assessment.

#### 1 Data Acquisition

I obtained a kidney scRNA-seq dataset from the HuBMAP Data Portal in .h5ad format.

The dataset contained 10,999 cells and 60,286 genes.

Since the dataset used Ensembl gene IDs, I needed to map them to gene symbols before annotation.

#### 2  Data Preprocessing in Python

Loaded the .h5ad file using anndata.

Inspected the structure of the dataset to confirm the presence of necessary data layers.

Set the correct expression matrix to .X and cleaned unnecessary metadata to avoid errors in downstream analysis.

Identified that the dataset contained Ensembl gene IDs instead of HGNC gene symbols, which required conversion.

#### 3 Gene Symbol Conversion

Used the MyGeneInfo API to map Ensembl Gene IDs to HGNC symbols.

Wrote a script to:

Query MyGeneInfo for each gene in adata.var_names.

Replace Ensembl IDs with corresponding HGNC symbols where available.

Retain Ensembl IDs for genes without mappings to avoid data loss.

Validated the conversion by printing the first few adata.var_names to confirm that gene symbols were correctly assigned.

#### 4 Normalization & Log Transformation

Normalized the dataset to 10,000 counts per cell using scanpy.pp.normalize_total().

Applied log1p transformation to stabilize variance in gene expression values.

Created a raw copy of the dataset for reference.

#### 5 Cell Type Annotation Using CellTypist

Downloaded the 'Immune_All_Low.pkl' pre-trained model from CellTypist.

Ran CellTypist annotation using:

prediction = celltypist.annotate(
    adata,
    model=model,
    majority_voting=True
)

Stored the predicted labels in the adata.obs['cell_type'] column.

Analyzed the output by printing the cell type distribution.

#### 6 Visualization of Cell Type Distribution

Created a bar plot showing the number of cells per predicted cell type.

Generated UMAP plots to visualize spatial separation of cell types.

#### 7 Saving and Exporting Results

Saved the final annotated .h5ad file:

adata.write("E:/annotated_expr_celltypist.h5ad")

Exported a CSV file containing the predicted labels for easier inspection:

adata.obs[['cell_type']].to_csv("E:/cell_annotations.csv")

🚧 Challenges & Issues Encountered

While working on this project, I encountered several challenges that required troubleshooting and problem-solving.
### 1 Missing Libraries & Dependencies

Faced multiple missing package errors (Seurat, Azimuth, SeuratDisk).

Installed them using CRAN and Bioconductor but ran into additional dependency issues.

#### 2 Azimuth Reference Issues

The initial plan was to annotate using Azimuth (R-based tool), but I could not find a suitable kidney reference file.

Attempted to install kidneyref from SeuratData but received 404 errors.

Also tried LoadReference("human_kidney.h5seurat"), but the reference file was unavailable.

As a result, I had to abandon Azimuth and pivot to CellTypist.

#### 3 Gene Symbol Compatibility Issues

The dataset used Ensembl gene IDs, but CellTypist requires HGNC symbols.

Implemented MyGeneInfo conversion but encountered missing mappings.

Addressed this by keeping unmapped genes in their original form to avoid data loss.

#### 4 Time Constraints

I found out about this assessment quite late because I was traveling to India for a wedding.

Due to limited time, I was unable to complete the R-based approach and focused on Python instead.

If I had more time, I would have implemented the Azimuth-based approach in R, given my prior experience with R and Seurat.

### Final Thoughts

While I would have loved to complete this assignment in R, the time constraints made it challenging. Given my prior experience in R, I am confident that if given more time, I could have executed the full workflow using Azimuth.

However, I adapted quickly and found a solution-oriented approach by using CellTypist, which allowed me to complete the annotation process successfully in Python.

I truly believe that my ability to learn fast, problem-solve, and adapt to unexpected challenges makes me a strong candidate for this role. I am eager to contribute and confident in my ability to perform the necessary duties efficiently.



