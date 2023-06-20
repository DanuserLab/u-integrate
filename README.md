# Related paper

- This repository is for initial distribution of the entire codes to generate the results in the manuscript, [**Granger-causal inference of the lamellipodial actin regulator hierarchy by live cell imaging without perturbation**](https://doi.org/10.1016/j.cels.2022.05.003), *Cell Systems*, 2022, 13(6):471-487.e8, written by Jungsik Noh, Tadamoto Isogai, Joseph Chi, Kushal Bhatt, [Gaudenz Danuser](https://www.danuserlab-utsw.org/).
- This initial distribution presents all the Matlab scripts for the pipeline, including a user-friendly Graphical-User-Interface to run the image processing of 2-channel live cell videos, cross-correlation analysis with molecular activities near the cell edge, and the Granger-causality inference pipeline. 

# Documents

- [Introduction to Granger-causality analysis of lamellipodia](/code/GrangerCausalityAnalysis/Doc/20210325_GC25min.pdf)
- [Introduction to Cross correlation Analysis](/code/mapDiagnosticsXcorrCurves/doc/readme_mapDDX.pdf)
- [Introduction to Fluctuation profiling around edge motion events](/code/FluctuationProfilingAroundEdgeMotionEvents/doc/readme_FPAEME.pdf)

# Software requirement

- The pipeline was built on Matlab R2020a.

# Running the GCA pipeline via scripts and example output

- A single Matlab script to run GrangerCausalityAnalysis (GCA) pipeline is [Pipeline_GCA_2chMovies_ch1ch2.m](code/GrangerCausalityAnalysis/Pipelines/Pipeline_GCA_2chMovies_ch1ch2.m). 
- Workflow
  1. (Part 1) Extract subcellular edge motion profiles and molecular activities from two channel live cell movies, using 'Windowing and Protrusion Package' (included in [./software](software)).
  2. (Part 2) Run [Pipeline_GCA_2chMovies_ch1ch2.m](code/GrangerCausalityAnalysis/Pipelines/Pipeline_GCA_2chMovies_ch1ch2.m) to implement Cross correlation analysis, Fluctuation profiling, and Granger-causality analysis. The sub-functions for these three modules are located under the [code](code/) folder. 
  3. Specific scripts for the above (Part 1) are described at [a manual pdf file](code/GrangerCausalityAnalysis/Doc/GCA_Pipeline_Manual_fromStartToFinalOutput.pdf).
  4. Example output of the pipeline is presented under the [example_output](example_output/) folder. For example, GC pathway diagram output for a toy dataset can be found [here](example_output/cropped-Tada200526_mDia1cr-mNG_Actin-Halo_downSampBilinear/ML_6by6ct2_plusSeg/MLmovies/GCA_3Variables_LF20fr_PL_GCA_2chMov_ewma0p5_Actin_mDia1cr_lL1wL1tL10). 

# Running the GCA pipeline via GUI

- The GUI can be simply installed by downloading the [./software](software) folder into your Matlab path.
- The following steps with GUI will generate inferred GC network diagrams from your input of live cell 2-channel videos.
- Suppose you have three videos of single cells where you imaged two protein activities, for example, actin and Arp2/3, as follows.

> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/inputFolders.png" /> 

- Launch the GUI using the following command in Matlab. Then click 'New' to create movieData objects for individual videos. 
```
>> movieSelectorGUI()
```
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/movieSelectorGUI_newMD.png" height="500"/> 

- Once managed to make movieData objects, click 'Save as movie list' in the movieSelectorGUI panel to create a movieList object.
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/makeML.png" height="500"/> 

- Now movieData objects and a movieList object are created. They will be the input of Segmentation, Protrusion, Windowing, Cross-correlation, Granger-causality analysis packages.
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/afterML.png" height="500"/> 

- Go to 'Windowing' package to run from segmentation to windowing. 
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/selectWindowing.png" height="500"/> 

- Run from Step 1: Generate Summation Channel to Step 8: Window Sampling for all the videos. In Step 2: Segmentation > MSA Segmentation process, check 'Use Output from Summation Channel', when the segmentation is implemented to the sum images as in the paper. 
- A user should make sure that the segmentation results are correct and accurate. You should adjust the segmentation parameters like 'tightness' and 'refinement radius' according to the imaging dataset after watching the segmentation results. 
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/MSApar_useOutputFromSummationChannel.png" height="500"/> 

- In Step 5: Protrusion process, set 'Mask process' to be 'Mask refinement' so that protrusion vectors are computed based on the refined masks.
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/protrusion_setrefinement.png" height="500"/> 

- In Step 6: Windowing process, set proper window sizes. In the paper, it was set to 6 pixels (720 nanometers).
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/setWindowSieze.png" height="500"/> 

- Once implemented the windowing package, go to 'Granger-Causality Analysis' package from the movieSelectorGUI panel.
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/selectGCApackage.png" height="500"/> 

- The first step of the GCA pacakge is to set up global parameters such as channel indexes, channel name, whether to apply low-frequency subtraction, number of layers to be analyzed, etc (See the paper).
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/GCApar.png" height="500"/> 

- You can run the rest steps of the GCA pipeline using the local parameters that are determined based on the global parameters set in Step 1. 
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/runGCA.png" height="500"/> 

- Output folders of the GCA pipeline
> <img src="/code/GrangerCausalityAnalysis/Doc/GUI_figures/GCAoutput.png" /> 



# Contact

Jungsik Noh (jungsik.noh@utsouthwestern.edu), Qiongjing (Jenny) Zou (Qiongjing.Zou@utsouthwestern.edu)

----------------------
[Danuser Lab Website](https://www.danuserlab-utsw.org/)

[Software Links](https://github.com/DanuserLab)
