% QZ change all 'lagMax0', 40 to 'lagMax0', 4
% QZ change all 'Arp3' & 'Apr3'(typo) to 'Actin' b/c wrong labelling.
% Delete all local functions.
% Qiongjing (Jenny) Zou, Aug 2018

% Updated on 11/20/2018
% add ML level calculation for LB Test and ActivityMapDescription.
% separate mapDDX to step 2 and 3.
% delete the Xcorr ch1 v.s. ch2 (chN v.s. ChN+1, when N>=1).

%% example_pipeline_DX_FPAEME_2chan_LB.m
%% The only input for this pipeline is the movieList object loaded to workspace.
% Jungsik Noh, 2018/05/08

load('/work/bioinformatics/s184919/Data/Jungsik/MLwith2Chan_premadeML_script_noMyproc/ML2mov/movieList.mat')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1 -- LB test: (same as 1chan_LB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Identifying quiescent windows by using Ljung-Box test

for i=1:numel(ML.movieDataFile_)
    %%  get movieData
    load(fullfile(ML.movieDataFile_{i}));
    pause(1)                                                
                                                    

    %%  specify outputDir

    outDirName = 'mapDescriptives_Vel_LB';
    figuresDir = fullfile(MD.outputDirectory_, outDirName);


    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end

    %%  run mapDescriptives_Vel_LB with options

    mapDescriptives_Vel_LB(MD, figuresDir, 'impute',1,'movingAvgSmoothing',1, ...
        'figFlag','off', 'omittedWindows', [], 'subFrames', subFr, 'topograph', 'on')

    close all
end

%% ML level
LBdirName = outDirName;
MLoutDirName = 'Summary_QuiescentWindow';
MLsummary_quiescentWindow(ML, LBdirName, 'ratioThreshold', 0.5, 'outDirName', MLoutDirName)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2 -- map Descriptive, (Diagnotics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  map Diagnostic and XCorrelation analysis (DX) for only non-quiescent windows

% for i=1:numel(ML.movieDataFile_)
% % movieList, MDindex, maxLayer, chname, index chan, mapDesc outDirName, XCorr outDirName
%     ML_DX_1chan_LB(ML, i, maxLayer, 'Actin', 1, 'mapDescriptives_LB', 'mapCrossCorr_LB')
%     ML_DX_1chan_LB(ML, i, maxLayer, 'mDia1', 2, 'mapDescriptives_LB', 'mapCrossCorr_LB')
%     ML_DX_xcorr2chan_LB(ML, i, maxLayer, 'mDia1', 'Actin', 2, 1, 'mapCrossCorr_LB')
% end

maxLayer = 2

ML.getMovies();
MDexp = ML.getMovie(1);
for j = 1:numel(MDexp.channels_)
    iChan = j;
    if j == 1 
        chNametag = 'Actin';
    elseif j == 2
        chNametag = 'mDia1';
    end

    for i=1:numel(ML.movieDataFile_) % MD Index
        %% get movieData at the MDindex-th location
        load(fullfile(ML.movieDataFile_{i}))
        pause(1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% At MD level,calculate omitWin and subFr for all following processes.
        %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
        %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
        velAnalName = 'mapDescriptives_Vel_LB';
        indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
        raw = load(indPath);
        
        omitWin = find(raw.indActive == 0);
        disp(omitWin)
        
        %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
        fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
        if ~isempty(fInfo)
            load(fullfile(MD.outputDirectory_, 'subFr.mat'));
            disp(['== length of subFr:', num2str(numel(subFr))])
        else
            subFr = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%  specify outputDir
        mapDescriptivesDirName = 'mapDescriptives_LB';
        figuresDir = fullfile(MD.outputDirectory_, mapDescriptivesDirName);
        
        %%  velocity Descriptives 
        % This part is exactly the same for chan1 and chan2; no need to run
        % twice in process.
        
        iChan0 = 0;
        chan0Name = 'Vel';
        chan0Title = 'Velocity (nm/sec)';
        
        mapDescriptives_OneChan(MD, iChan0, 1, chan0Name, chan0Title, figuresDir, ...
            'adf',1,'impute',0, 'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 100, 'topograph', 'on')
        
        
        %%  Chan Descriptives
        chanName = chNametag;
        chanTitle = chanName;
        
        mapDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, figuresDir, ...
            'adf', 1, 'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr,'numPerm', 100, 'topograph', 'on', 'WithN', 0)
        
        close all
    end
end

%% ML level
mapDescDirName = mapDescriptivesDirName;
for j = 1:numel(MDexp.channels_)
    iChan = j;
    if j == 1 
        chNametag = 'Actin';
    elseif j == 2
        chNametag = 'mDia1';
    end
    MLoutDirName = 'Summary_ADFtest';
    MLsummary_ADFtest(ML, mapDescDirName, iChan, chNametag, 'outDirName', MLoutDirName)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 3 -- Cross correlation analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:numel(MDexp.channels_)
    iChan = j;
    if j == 1 
        chNametag = 'Actin';
    elseif j == 2
        chNametag = 'mDia1';
    end

    for i=1:numel(ML.movieDataFile_) % MD Index
        %% get movieData at the MDindex-th location
        load(fullfile(ML.movieDataFile_{i}))
        pause(1)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% At MD level,calculate omitWin and subFr for all following processes.
        %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
        %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
        velAnalName = 'mapDescriptives_Vel_LB';
        indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
        raw = load(indPath);
        
        omitWin = find(raw.indActive == 0);
        disp(omitWin)
        
        %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
        fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
        if ~isempty(fInfo)
            load(fullfile(MD.outputDirectory_, 'subFr.mat'));
            disp(['== length of subFr:', num2str(numel(subFr))])
        else
            subFr = [];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%  specify outputDir of correlation analysis
        mapXCorrDirName = 'mapCrossCorr_LB';
        figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);
        
        %% xcorr chan vs velocity
        mapXcorrCurvePermutation_Vel(MD, iChan, chNametag, maxLayer, figuresDir, ...
            'impute', 0, 'movingAvgSmoothing',0, 'omittedWindows', omitWin, 'subFrames', subFr, 'WithN', 0, 'numPerm', 100, 'lagMax', 4)
        
        close all
        
    end

end

%% xcorr ch1, ch2 at MD level
% for i=1:numel(ML.movieDataFile_)
%     % get movieData at the MDindex-th location
%     load(fullfile(ML.movieDataFile_{i}))      %ML.getMovies();;
%     pause(1)                                                %MD = ML.getMovie(MDindex);
% 
%     disp('===========================')
%     disp(['======= MDindex = ', num2str(i), ' ======='])
%     disp('===========================')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% At MD level,calculate omitWin and subFr for all following processes.
%     %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
%     %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
%     velAnalName = 'mapDescriptives_Vel_LB';
%     indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
%     raw = load(indPath);
% 
%     omitWin = find(raw.indActive == 0);
%     disp(omitWin)
% 
%     %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
%     fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
%     if ~isempty(fInfo)
%         load(fullfile(MD.outputDirectory_, 'subFr.mat'));
%         disp(['== length of subFr:', num2str(numel(subFr))])
%     else
%         subFr = [];
%     end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %%  outputDir
%     mapXCorrDirName = 'mapCrossCorr_LB';  % same folder as xcorr chan vs velocity
%     figuresDir = fullfile(MD.outputDirectory_, mapXCorrDirName);
% 
% 
%     %% xcorr ch1, ch2  %% NOTE iChan1 is 2, and iChan2 is 1 !
%     iChan1 = 2;
%     iChan2 = 1;
%     chNametag1 = 'mDia1';
%     chNametag2 = 'Actin';
%     chan1Name = chNametag1;
%     chan2Name = chNametag2;
% 
%     mapXcorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name, maxLayer, figuresDir, ...
%        'impute', 0, 'parpoolNum', 4, 'WithN', 0, 'numPerm', 100, 'omittedWindows', omitWin, ...
%        'movingAvgSmoothing', 0, 'subFrames', subFr) %QZ I changed 'parpoolNum' from 12 to 4.
%         
%     close all
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Xcorr analysis summary at the ML level

maxLayer = 2
% NOTE different out put file!
MLsummary_XcorrCurvesVelAcf(ML,1,0,'Actin','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 4, ...
    'outDirName', 'Xcf_ch1ch0_LB_lagMax4')
close all
%
MLsummary_XcorrCurvesVelAcf(ML,2,0,'mDia1','Vel', maxLayer, ...
    'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 4, ...
    'outDirName', 'Xcf_ch2ch0_LB_lagMax4')
close all
%
% MLsummary_XcorrCurvesVelAcf(ML,2,1,'mDia1','Actin', maxLayer, ...
%     'mapDescriptives_LB', 'mapCrossCorr_LB', 'lagMax0', 4, ...
%     'outDirName', 'Xcf_ch2ch1_LB_lagMax4')
% 
% close all
%%  end: DX


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 4 -- FPEAME:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  Fluctuation Profiling Around Edge Motion Events (FPAEME)
% this section has no interaction btw chan1 and chan2, 
% in MLsummary_FluctuationProfiling there is some interaction/calculation btw MD1 and MD2.

maxLayer = 2

   %%%%%%%%%%%
    % iChan = 1; % NOTE only diff btw channels
    % chNametag = 'Actin'; % NOTE only diff btw channels

ML.getMovies();
MDexp = ML.getMovie(1);
for j = 1:numel(MDexp.channels_)
    iChan = j;
    if j == 1
        chNametag = 'Actin';
    elseif j == 2
        chNametag = 'mDia1';
    end

    %%  ML_phaseDescriptives
    tic


    for i=1:numel(ML.movieDataFile_)
        % check smParam and minimumRunLength within ML_ function.

        %% get movieData at the MDindex-th location
        load(fullfile(ML.movieDataFile_{i}))      
        pause(1)    

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% At MD level,calculate omitWin and subFr for all following processes.
    %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);

    omitWin = find(raw.indActive == 0);
    disp(omitWin)

    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%  specify outputDir
        outDirName = 'phaseDescriptives';
        figuresDir = fullfile(MD.outputDirectory_, outDirName);

        %% phaseMasking to find a proper smoothing parameter

        smParamTh = 0.4

        [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, 'minimumRunLength', 5, ...
                            'subFrames', subFr);

        %% specify sampling bandwidth
        samplingBw = 20

        %%  Retraction 
            Mask0 = retMask1;

            %%  Chan Descriptives 
            
            chanName = [chNametag, '-ret'];
            chanTitle = chanName;

            %
            phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
                    'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)

        %%  Protrusion
            Mask0 = protMask1;

            %%  Chan Descriptives

            chanName = [chNametag, '-prot'];
            chanTitle = chanName;

            %
            phaseDescriptives_OneChan(MD, iChan, maxLayer, chanName, chanTitle, Mask0, samplingBw, figuresDir, ...
                    'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)


        %%  phase MaxMinVel 

        chanName = chNametag;
        chanTitle = chanName;

        phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, chanTitle, ...
                smParamTh, samplingBw, figuresDir, ...
                'WithN', 0, 'impute', 1, 'omittedWindows', omitWin, 'subFrames', subFr, 'movingAvgSmoothing', 0)

        close all

    end

    toc
    %  FPAME summary at the ML level
    outDirName = ['Summary_PhaseDesc_' chNametag '_Chan', num2str(iChan)];

    MLsummary_FluctuationProfiling(ML, chNametag, maxLayer, 'phaseDescriptives', ...
        outDirName, 'lagMax0', 5)
    
    close all

end


% %  ML_phaseDescriptives
% tic

% for i=1:numel(ML.movieDataFile_)
%     % check smParam and minimumRunLength within ML_ function.
%     ML_phaseDescriptives_LB(ML, i, maxLayer, 'mDia1', 2, 'phaseDescriptives0p5mRL5_LB')
% end
% toc
% %  FPAME summary at the ML level

% outDirName = ['Summary_PhaseDesc_' chNametag];
% MLsummary_FluctuationProfiling(ML, 'mDia1', maxLayer, 'phaseDescriptives', ...
%     'phaseDesc_mDia1_LB_0p5mRL5_lagMax20', 'lagMax0', 20)



%%  end: FPAME


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% (not on GUI) -- alongDepth:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%  ML_meanCV_alongDepth_LB
% here only calculated for chan1 'Actin', can do for chan2 'mDia1', too.
maxLayer = 5;

ML.getMovies();  % Get the movies from a movie list
MDexp = ML.getMovie(1);  % Get the movies from a movie list
for j = 1:numel(MDexp.channels_)
    iChan = j;
    if j == 1 
        chNametag = 'Actin';
    elseif j == 2
        chNametag = 'mDia1';
    end

    tic
    for i=1:numel(ML.movieDataFile_)

        %% get movieData at the MDindex-th location
        load(fullfile(ML.movieDataFile_{i}))      %ML.getMovies();;
        pause(1)                                                %MD = ML.getMovie(MDindex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% At MD level,calculate omitWin and subFr for all following processes.
    %%  exclude quiescent windows which are classified by mapDescriptives_Vel_LB() in advance.
    %%  outDir of mapDescriptives_Vel_LB() should be the inputDir (=velAnalName) here.
    velAnalName = 'mapDescriptives_Vel_LB';
    indPath = fullfile(MD.outputDirectory_, velAnalName, 'indActive_windowIndex.mat');
    raw = load(indPath);

    omitWin = find(raw.indActive == 0);
    disp(omitWin)

    %%  Only if MD.outputDirectory_ contains 'subFr.mat', analysis is done for the specified subframes.
    fInfo = dir(fullfile(MD.outputDirectory_, 'subFr.mat'));
    if ~isempty(fInfo)
        load(fullfile(MD.outputDirectory_, 'subFr.mat'));
        disp(['== length of subFr:', num2str(numel(subFr))])
    else
        subFr = [];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%  outputDir
        outDirName = 'mapDescriptives_alongDepth';
        figuresDir = fullfile(MD.outputDirectory_, outDirName);

        %%  
        chanName = chNametag; 

        mapDescriptives_meanCV_alongDepth(MD, iChan,maxLayer, chanName,figuresDir, ...
            'impute',0,'WithN', 0, 'subFrames', subFr, 'omittedWindows', omitWin)

        close all

    end
    toc

    %  MLsummary_alongDepth
    analFolderName = outDirName;
    outDir = ['topolayers_', 'Chan', num2str(iChan)];
    MLsummary_alongDepth(ML, iChan, maxLayer, chNametag, analFolderName, ...
         'outDirName', outDir)
    
    close all
    
end
%%  EOF
