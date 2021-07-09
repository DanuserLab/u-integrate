% A sample script to run Cross Correlation and Fluctuation Profile Package GUI by script.
% ML is a movielist with 2 movieData and 6 processes in windowing package already run.
% Qiongjing (Jenny) Zou, Oct 2018

load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_allScript/ML2mov/movieList2chan.mat')
ML.reset;


package = XcorrFluctuationPackage(ML);
ML.addPackage(package);

ML.getPackage(1).createDefaultProcess(1); 
%     funParams = ML.getProcess(1).funParams_;
%     funParams.figFlag = 'on';
%     ML.getProcess(1).setParameters(funParams);
ML.getPackage(1).getProcess(1).run();

ML.getPackage(1).createDefaultProcess(2); % default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(2).funParams_;
    funParams.omittedWindows = true;
%     funParams.parpoolNum = 6;
%     funParams.chanName = {'Actin' 'mDia1'}; % rename channel names.
    ML.getProcess(2).setParameters(funParams);   
ML.getPackage(1).getProcess(2).run(); 

close all

ML.getPackage(1).createDefaultProcess(3); % default does not use LB Test result.
    % run below 3 lines to use LB Test result:
    funParams = ML.getProcess(3).funParams_;
    funParams.omittedWindows = true;
%     funParams.chanName = {'Actin' 'mDia1'}; % rename channel names.
    ML.getProcess(3).setParameters(funParams);  
ML.getPackage(1).getProcess(3).run();
   
close all