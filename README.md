# Yellow-birch-photosynthesis

This project directory serves as basis for the Yellow-Birch-Photosynthesis project. It was created on September 7th 2023 by Lukas Van Riel. What follows is a project description, a description of the different folders in this root folder and the used naming convention. More details can be found on the OSF project page: TODO INSERT LINK.

### Project description

This is a field study that seeks to examine the present and future effects of climate change on the physiology of temperate tree species near their northern range limit in Quebec, Canada. Yellow birch (Betula alleghaniensis Britton) will function as exemplary species. This research aims to answer the following questions in particular:

-   How does the photosynthetic capacity of yellow birch saplings and adults vary over a latitudinal gradient?
-   How does the fraction of carbohydrates allocated for growth vs. reproduction vary with latitude for yellow birch?

The study area consists of five evenly spaced sites, each consisting of three plots, along a latitudinal gradient between 46 °N and 48 °N located close to 74.5 °W. All plots will be circular with an area of 400 m2 and are located within 25 km of their respective site centroid. At each plot, gas exchange measurements of five adult and five sapling (diameter at breast height = 1-9 cm) trees will be taken with the LI-6800 Portable Photosynthesis System (Li-Cor Inc., Lincoln, NE, USA) using the rapid A/Ci curve (RACiR; Stinziano et al., 2017) method. These measurements will then be converted to relevant photosynthesis parameters using physiological models in order to see whether climate is currently the primary determinant of the species' northern limit as well as project the performance under future climatic conditions. In addition, tree core samples will be taken to determine the radial growth rates of the individuals as well as the spectral properties of each sampled leaf with a field spectroradiometer (PSR- 3500 Spectral Evolution Inc.) to complement the RACiR measurements. The trees of each plot will be mapped, the plant community will be recorded and the plot canopy openness measured using a LAI- 2200C Plant Canopy Analyzer (Li-Cor Inc., Lincoln, NE, USA). A soil sample will be taken that will be analysed in the lab to consider any soil variability. The other abiotic variables such as ambient temperature and vapour pressure deficit will be recorded by the LI-6800 system during measurements.

### Directory structure

00_rawdata: contains all the raw data that was collected. 

01_Scripts: contains all the necessary scripts to analyse the data.

02_outdata: contains the output data used to create figures, results...

03_figs: contains all the figures that were created by the analysis.

04_manuscripts: contains the versions of the preregistration and final manuscript

### Naming convention

00_rawdata: LastName_MeasurementSystemName_YYYYMMDD.csv

01_Scripts: LastName_Scriptname.R

02_outdata: LastName_YYYYMMDD_TableName.csv

03_figs: LastName_Parameter/Type_FigureType.png

04_manuscripts: LastName_FileType.type


<br>


*Last update: 2023-09-26 by Lukas Van Riel.*
