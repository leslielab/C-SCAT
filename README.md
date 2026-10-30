This repository contains examples of the MATLAB C-SCAT code used in "Quantitative Label-free Characterization of Nanoparticles and Subcellular Structures in Cells Using Confocal Scattering Microscopy", soon to be published in npj Imaging. The current example codes and data includes particle on glass, particle in cells, and structures in cells. Note that that the amount of example data per example case is limited compared to the what is shown in the manuscript due to the size of the microscopy videos.

Requirements: The code is based on MATLAB 2022 with 'Image processing toolbox' and 'Curve fitting toolbox'.

The code consists of two main parts, where one is for particle analysis and the other is for cell segment analysis. For the particle analysis, use the code 'Main_particleAnalysis.m', which consists of an analysis part and an plot part. For the cell segment analysis, use the code 'Main_CellSegment.m', which goes through the process of selecting and analysis cellular segments.

To adapt the code so that it can run with your data, you first need add a way to loading your data into MATLAB (in the example the data is .mat files). In addition to this, the image processing is currently based on that the data has 4 simultanous imaging channels, where that can be changed depending on the your data. There are also some thresholds for particle detection etc., where those details can be found when you open the '.m' files.

Since '/' and '' are different when using Mac or Windows, there is an option in the code to switch between the two different ways of loading the data.

The overall data pipelines for particle and cell structure analysis are summarized in the figures below. For more details, see Supporting Information of "Quantitative Label-free Characterization of Nanoparticles and Subcellular Structures in Cells Using Confocal Scattering Microscopy".
