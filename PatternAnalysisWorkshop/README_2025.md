# Tutorial: Evaluating Digital Soil Maps by their patterns -- Workshop Bangalore 2025

This *tutorial*, written as a Quarto document, presents methods to evaluate the spatial patterns of the spatial distribution of soil properties and map units as shown in gridded maps produced by digital soil mapping (DSM). It compares patterns from different DSM products to each other, and to spatial patterns known from detailed field surveys or known to local experts but not represented (yet) on maps. Methods include whole-map statistics, visually identifiable landscape features, level of detail, range and strength of spatial autocorrelation, landscape metrics (Shannon diversity and evenness, shape, aggregation, mean fractal dimension, and co-occurrence vectors), and spatial patterns of property maps classified by histogram equalization or user-defined cutpoints. The tutorial also shows how to use patterns within a window to partition a soil landscape into zones with similar patterns. This workshop uses an example window from SoilGrids v2.0 , but the methods are applicable to any gridded DSM product.

Quarto source (load into R Markdown): `PatternAnalysisWorkshopTutorial_DSM2025.qmd`

Supporting `BibTeX` bibliography file: `patterns.bib`

Compiled documents: `.pdf` and `.html` (with supporting folders)

R Markdown source to install required packages: `InstallWorkshopPackages_2025.R` 

R Script to download soil properties within a window from SoilGrids v2.0.

R Markdown sources to build a GeoTIFF from SoilGrids digital soil maps: `SoilGrids250_MakeRasterStack.Rmd`

