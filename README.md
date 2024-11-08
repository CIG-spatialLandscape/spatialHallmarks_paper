
<h1>
 <img align="right" src="map.png" alt="Example Image" width="150" height="150">
 Unraveling the spatial architecture<br />
  of Cancer Hallmarks
</h1>

#### Authors: Mustafa Sibai* & Sergi Cervilla*
\* Equal contribution

This is the updated code corresponding to our latest submitted manuscript (under-review), which is used for the collection of Cancer Hallmark gene sets and the computation of their activities in enhanced Visium datasets from 63 samples across 10 cancer types with all downstream analysis. Our previous version of this study correspond to this preprint: doi: https://doi.org/10.1101/2022.06.18.496114  <br />

## Outline
- SpaceRanger
- Processing
- Enhancing
- Hallmark gene collection
- Hallmark activity computation
- ESTIMATE computation
- Pan-Cancer heatmap
- CNV estimation
- Random Forest: Radar computation and RF models
- Plots
- Supplementary
- Utils (Annotation, Plot functions, compute real coordinates...) 

## Code purpose
- **Preprocess.R**: Preprocess samples, transform to SingleCellExperiment object and run clustering
- **EnhanceST.R**: Use BayesSpace enhancing algorithm in Visium Spatial transcriptomics data
- **GeneImputation.R**: Impute genes at sub-spot resolution
- **GeneCollection.R**: Obtain hallmark gene-lists using Pathway Commons
- **HallmarkScores.R**: Script to create a Seurat object from an enhanced SingleCellExperiment object and run AddmoduleScore on sub-spot resolution
- **mat2GCT.sh**: Tranform expression matrix to GCT format
- **ComputeESTIMATE.R**: Extract enhanced expression matrix and compute ESTIMATE scores at sub-spot resolution
- **SpAutocorrelation**: Compute Moran's for hallmark activities at sub-spot resolution
- **CNVexperiment.R**: Generate CNV clusters, run CNV experiment at spot level and transfer to sub-spot level
- **CancerRadar.R**: Compute TME Radar scores for Cancer sub-spots
- **TMERadar.R**: Compute Cancer Radar scores for TME sub-spots
- **RFCancer.R**: Generate Random Forest model to predict a given Cancer Hallmark in a given sample at sub-spot resolution
- **RFTME.R**: Generate Random Forest model to predict a given TME Hallmark in a given sample at sub-spot resolution
- **Utils** 
  - **CoordinatesEnhanced.R**: Function to compute the real distance between sub-spots
  - **SamplesMetadata.R**: Variables with metadata and full names for hallmarks and samples
  - **PlottingMod.R**: Modifications of BayesSpace and Seurat plotting to combine and plot high resolution sub-spot plots
  - **annotLook.txt**: Gene type annotation look up table
  - **header.txt**: Header template for GCT files
  - **sample_ID.txt**: table to match samples ID in circos plots
- **Plots** 
  - **CreateHiresFiles.R**: Create files to plot high resolution images with enhanced BayesSpace sub-spots
  - **PanCancerHeatmap.R**: Plot a heatmap representing and the hallmark activity within each ESTIMATE cluster
 
