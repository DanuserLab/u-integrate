# Related paper

- This repository is for initial distribution of the entire codes to generate the results in the manuscript, "[Inference of Granger-causal relations in molecular systems — a case study of the functional hierarchy among actin regulators in lamellipodia](https://www.biorxiv.org/content/10.1101/2021.05.21.445144v2)." at BioRxiv.
- This initial distribution presents all the Matlab scripts for the pipeline, which is not yet optimized for easy-to-use. We will soon deploy a user-friendly Graphical-User-Interface to run the GCA pipeline. 


# Accessing the codes and example output

- A single Matlab script to run GrangerCausalityAnalysis (GCA) pipeline is [Pipeline_GCA_2chMovies_ch1ch2.m](code/GrangerCausalityAnalysis/Pipelines/Pipeline_GCA_2chMovies_ch1ch2.m). 
- Workflow
  1. (Part 1) Extract subcellular edge motion profiles and molecular activities from two channel live cell movies, using 'Windowing and Protrusion Package' available at [a different repository](https://github.com/DanuserLab/Windowing-Protrusion).
  2. (Part 2) Run [Pipeline_GCA_2chMovies_ch1ch2.m](code/GrangerCausalityAnalysis/Pipelines/Pipeline_GCA_2chMovies_ch1ch2.m) to implement Cross correlation analysis, Fluctuation profiling, and Granger-causality analysis. The sub-functions for these three modules are located under the [code](code/) folder. 
  3. Specific scripts for the above (Part 1) are described at [a manual pdf file](code/GrangerCausalityAnalysis/Doc/GCA_Pipeline_Manual_fromStartToFinalOutput.pdf).
  4. Example output of the pipeline is presented under the [example_output](example_output/) folder. For example, GC pathway diagram output for a toy dataset can be found [here](example_output/cropped-Tada200526_mDia1cr-mNG_Actin-Halo_downSampBilinear/ML_6by6ct2_plusSeg/MLmovies/GCA_3Variables_LF20fr_PL_GCA_2chMov_ewma0p5_Actin_mDia1cr_lL1wL1tL10). 


# Documents

- [Introduction to Granger-causality analysis of lamellipodia](/code/GrangerCausalityAnalysis/Doc/20210325_GC25min.pdf)
- [Introduction to Cross correlation Analysis](/code/mapDiagnosticsXcorrCurves/doc/readme_mapDDX.pdf)
- [Introduction to Fluctuation profiling around edge motion events](/code/FluctuationProfilingAroundEdgeMotionEvents/doc/readme_FPAEME.pdf)


# Software requirement

- The pipeline was built on Matlab R2020a.


# Contact

Jungsik Noh (jungsik.noh@utsouthwestern.edu), Qiongjing (Jenny) Zou (Qiongjing.Zou@utsouthwestern.edu)
