

Manuscript tittle and authors

  The code of this folders it has been submitted joint with the manuscript entilted "Approximate Bayesian inference for multivariate point pattern analysis in disease mapping" to the Biometrical Journal. 

  The authors of this work are Francisco Palmí Perales (Francisco.Palmi@uclm.es, Corresponding author) and Virgilio Gómez Rubio (Virgilio.Gomez@uclm.es) from the Universidad de Castilla-La Mancha); Gonzalo López-Abente Ortega, Rebeca Ramis Prieto and Pablo Ferández Navarro from Instituto de Salud Carlos III and Jose Miguel Sánz AnquelaHospital from the Hospital Universitario Príncipe de Asturias.

  The code has been mainly developed by Francisco Palmí Perales and Virgilio Gómez Rubio. Should you have any question do not hesitate to send an email to any of this two email adressess Francisco.Palmi@uclm.es/ Virgilio.Gomez@uclm.es.


How to run the code

  R software and some of their packages have been used to develope this code. Thus you will need to install this sofware and the appropriate packages which you can check in the following lines. After loading all the packages needed the outpout of the sessionInfo() is the following: 

  R version 3.6.0 (2019-04-26)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

Random number generation:
 RNG:     Mersenne-Twister 
 Normal:  Inversion 
 Sample:  Rounding 
 
locale:
[1] LC_COLLATE=Spanish_Spain.1252  LC_CTYPE=Spanish_Spain.1252    LC_MONETARY=Spanish_Spain.1252 LC_NUMERIC=C                  
[5] LC_TIME=Spanish_Spain.1252    

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3         lattice_0.20-38       mapproj_1.2.7         maps_3.3.0            OpenStreetMap_0.3.4  
 [6] tmaptools_3.0         rgdal_1.4-8           ggmap_3.0.0           ggplot2_3.3.0         RColorBrewer_1.1-2   
[11] rgeos_0.5-2           spatstat.utils_1.17-0 spatstat_1.63-3       rpart_4.1-15          nlme_3.1-139         
[16] spatstat.data_1.4-3   maptools_0.9-9        INLA_20.03.17         foreach_1.5.0         Matrix_1.2-17        
[21] sp_1.4-1             

loaded via a namespace (and not attached):
 [1] httr_1.4.1          tidyr_1.1.0         viridisLite_0.3.0   splines_3.6.0       assertthat_0.2.1    pillar_1.4.4       
 [7] glue_1.4.0          digest_0.6.25       polyclip_1.10-0     colorspace_1.4-1    plyr_1.8.6          XML_3.99-0.3       
[13] pkgconfig_2.0.3     raster_3.1-5        purrr_0.3.4         scales_1.1.0        tensor_1.5          jpeg_0.1-8.1       
[19] tibble_3.0.1        mgcv_1.8-28         farver_2.0.3        ellipsis_0.3.0      withr_2.2.0         magrittr_1.5       
[25] crayon_1.3.4        deldir_0.1-25       lwgeom_0.2-3        foreign_0.8-71      class_7.3-15        tools_3.6.0        
[31] RgoogleMaps_1.4.5.3 lifecycle_0.2.0     stringr_1.4.0       munsell_0.5.0       Deriv_4.0           compiler_3.6.0     
[37] e1071_1.7-3         rlang_0.4.6         classInt_0.4-3      units_0.6-6         grid_3.6.0          dichromat_2.0-0    
[43] iterators_1.0.12    rstudioapi_0.11     rjson_0.2.20        goftest_1.2-2       labeling_0.3        bitops_1.0-6       
[49] gtable_0.3.0        codetools_0.2-16    abind_1.4-5         DBI_1.1.0           curl_4.3            R6_2.4.1           
[55] dplyr_0.8.5         KernSmooth_2.23-15  rJava_0.9-12        stringi_1.4.6       Rcpp_1.0.4.6        vctrs_0.3.0        
[61] sf_0.9-3            png_0.1-7           tidyselect_1.1.0   


Structure and the execution order of the files

  Apart from this README.txt, there are two main folders: case_study and simulation. The first folder contains the code used to obtain the results of the sections 5.1 to 5.4 of the manuscript and the second contains the code used to obtain the results shown in section 5.5. 

  The case_study folder contains sseveral numbered R files which will allow the user to fully reproduce the procedure of the manuscript using simulated data (the original data is not provided due to confidentiality reasons). Also, another R file can be found which contains the code to repoduce the figures of the manuscript (see below). Furthermore, the user will find two folders: 'Data' that contains the simulated data to test the code and 'Results' which will store the outcomes of the files. Finally, a file entitled "Non-publishable-folder.txt" which has been placed instead of the folder with the original data can be also found. 

  The 'simulation' folder contains several files. The simulated data is stored in sevral folders named SimulatedDataphi?.RData which have been produced using the code in the file Simulating-Data.R. Then, we can find four .R files (Running-M0-M1-SimulationStudy.R, Running-M2-M3-LinEff-SimulationStudy.R, Running-M2-M3-Rw1-SimulationStudy.R and Running-M2-M3-SPDE1D-SimulationStudy.R) which can be used to run the different models using the data of the simulation study. Also, a file with the necessary code to plot Figure 3 of the manuscript can be found (Code-Figure3.R). Finally, three more files can be found (Data-with-VariabilityAdded-1Controls-Less.txt, focos_aire_contam.txt and focos_metales_pesados.txt), the first one has been placed instead of the original data which can not be plublished. The other to correspond to the data of the location of the polluting industries considered in the study. 

Data source

  The original data cannot be shared due to confidentiality reasons (as it records the exact localtion of the cases), therefore it can not be distributed. However, we have provided a simulated data which is really similar to the original one in order to been able to reproduce the analysis using the code. This simulated data has been obtained from the original data by estimating the intensity of the original one and them simulating from this estamted intensity. 

Tables

  No code has been used to generate the tables which have been built by picking the results from the R output and creating a table with LaTeX. Table 1 of the manuscript is a summary of the different models. Table 2 sums up the DIC values of the simulation study results which have been obtained using the R files of the simulation folder. Table 3 is a summary of the results of the different models fitted which values have been extracted of the results obtained with the R files of the case_study folder. Finally, Table 4 is a summary the values of the model selection criteria which has been obtained from the results, specifically, inla objects can store these values where is specified in the code which can be seen in our code. 

Figures

  The code for reproducing the figures of the manuscript is in the file Figures-code which is located in the case_study folder. The code for ploting figure 3 can be found in the file Code-Figure3 which is located in the simulation folder. 


  
