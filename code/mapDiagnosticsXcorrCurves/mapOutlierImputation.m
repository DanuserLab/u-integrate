function [fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap, ...
    actmap_outlSc, LFnormalizingMap, feigCurve, VNout] ...
    = mapOutlierImputation(MD, iChan, maxLayer, varargin)
% mapOutlierImputation Get the input of MD, iChan and maxLayer and Give the
% corresponding 3-dim (window, layer, time frame) activity array in raw,
% outlier-removed (>5*sigma), and imputed forms.
%
% Usage:
%       [fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap]
%       = mapOutlierImputation(MD, iChan, maxLayer, varargin)
%
% Input:
%       MD          - a movieData object
%       iChan       - a channel index. If 0, the edge velocity is analyzed.
%                     iChan = 1x indicates the differenced map (X_t =
%                     X_{t-1}) of channel x.
%                     iChan = 2x indicates movmedian normalized (LFnormalize).
%                     iChan = 3x indicates LowFreq normalizing Map.
%                     iChan = 1xx indicates common factor-based normlized.
%                     iChan = 2xx indicates the CommonFactor normalizing
%                     Map.
%       maxLayer    - maximum layer to be analyzed
%
% Output:
%       fname0      - a suggested file name for the map. ['Chan', num2str(iChan)]
%       MDpixelSize_    - pixelSize of the movieData
%       MDtimeInterval_ - time interval of the movieData
%       wmax        - number of windows per layer or nSliceMax_
%       tmax        - number of time frames or MD.nFrames_
%       rawActmap   - 3-dim'l array (window, layer, time frame) of raw
%                   activities.
%       actmap_outl - 3-dim'l array of outlier-removed activities. Outliers
%                   are determined for each layer based on (|Z-score| > 5).
%       imActmap    - 3-dim'l array of activities in which missing values
%                   are imputed by using k-nearest neighbor method. If
%                   'impute' option is false, imActmap is the same as
%                   actmap_outl.
%       LFnormalizingMap     - if LowFrequency normalization with movMedian is used, then the baseline
%                          activity used in the normaliziation. Otherwise NaN.
%
% Option:
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function. Default is true.
%       WithN       - if true, it uses an alternative windowSampling result
%                   which is obtained by sampleMovieWindowsWithN.m and includes number
%                   of pixels for each windows. Default is false.
%       omittedWindows
%                   - window index in which activities will be replaced by
%                   NaN. Default is null.
%       subFrames
%                   - specified frames will be only used.
%       movingAvgSmoothing
%                   - input activities are minimally convolved with just
%                   adjacent activities, that is, by using 3 by 3 patch.
%
% Updated:
% J. Noh, 2019/05/31. Add a 'VNtlagMax' option for Vel-Normalization.
% J Noh. 2019/05/16. CF normalization is implemented before smoothing.
% J Noh. 2019/05/11. Introduce 'figuresDir', an option for the output directory
% to deal with complex preprocessing output. It now saves map output as an .mat
% file. Now 'CommonFactorNormAddCh' maps are loaded from 'figuresDir'.
% J Noh, 2019/05/09. CFnormMap is now in the Z-score scale for better
% visualization.
% J Noh, 2019/05/04. 'EWMA' specifies lambda value for exponentially
% weighted moving average smoothing. Default is 1 (raw data). It is also
% applied to CFnormalization.
% J Noh, 2019/04/30. Add 'ratio' chan (4xx) and a 'baseOfRatio' option.
% J Noh, 2019/02/28. Add an option, 'factoranMethod' (internally used before).
% J Noh, 2018/11/09. 'numConsecutiveFrames' -> 'factoranMethod'.
% J Noh, 2018/10/30. Option 'numConsecutiveFrames' for CommonFactor analysis.
% J Noh, 2018/10/28. Move 'LFnormalize', 'CFnormalize' options inside of
% mapOutlierImputation().
% J Noh, 2018/10/26. Movmedian normalized activity has now % unit (.*100).
% p.normalize => p.LFnormalize. Add p.LFof and p.CFof.
% J Noh, 2018/05/30. Add Common Factor (CF) analysis-based normalization.
% J Noh, 2018/05/28. Add 'outlSigma' option.
% J. Noh, 2018/05/27. Normalization using movmedian.
% Jungsik Noh, 2018/02/05.
% J Noh, 2017/10/11, raw activities can be smoothed. New option is
% 'movingAvgSmoothing'.
% J Noh, 2017/08/26, To deal with differenced channels.
% Jungsik Noh, 2017/05/23
% Jungsik Noh, 2016/10/24


ip = inputParser;
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('movingAvgSmoothing', false);
%ip.addParameter('movMedFrameSize', 30);
ip.addParameter('outlSigma', 5);
%ip.addParameter('LFnormalize', false);
ip.addParameter('movMedFrameSize', nan);
%ip.addParameter('CFnormalize', false);    %
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('factoranMethod', 33);  %%%% 33: 3 frs, 4 factors, 32: 3frs, 2 factors
ip.addParameter('baseOfRatio', nan);  % base chan for ratio comp. default is ch1.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing.
ip.addParameter('figuresDir', []);  % for output
ip.addParameter('chanName', []);    % necessary for output
ip.addParameter('VNtlagMax', 10);   % maximum lag (frames) for VelNorm with SPAR


parse(ip, varargin{:});
p = ip.Results;

%%  p.LFxxx according to iChan

p.Ratio = 0;
p.CFnormalize = 0;

if (iChan >= 400) && (iChan < 500); p.Ratio = 1;
end

if (iChan >= 100) && (iChan < 300)
    p.CFnormalize = 1;
end
% 522 is e.g. CFN.R.LFN.Arp3
if (iChan >= 500) && (iChan < 600); p.Ratio = 1; p.CFnormalize = 1;
end

if (iChan >= 200) && (iChan < 300)
    p.CFof = 1;
else
    p.CFof = 0;
end

%% p.VelNorm
p.VelNorm = 0;
% eg, VN.LFN.Actin
if (iChan >= 700) && (iChan < 800); p.VelNorm = 1;
end

%%
iChan_tmp = mod(iChan, 100);

if (iChan_tmp >= 20) && (iChan_tmp <40)
    p.LFnormalize = 1;
else
    p.LFnormalize = 0;
end


%%  getting Maps from channels (0 indicates the edge velocity)

disp(['==== Channel index: ', num2str(iChan), ' ===='])


%iChan_tmp = iChan;
%if p.CFnormalize
%    iChan_tmp = mod(iChan, 100);
%end

% 2017/08 p.derivative
% 2018/05 p.derivative, p.movMedNormalize modified
% 2018/10 p.LFnormalize is modified

if (iChan_tmp >= 10) && (iChan_tmp < 20)
    p.derivative = 1;
else
    p.derivative = 0;
end

mdChan = mod(iChan_tmp, 10);

if (iChan_tmp >= 30) && (iChan_tmp < 40)
    p.LFof = 1;
else
    p.LFof = 0;
end


%if (iChan_tmp >= 10) && (iChan_tmp < 20)
%    p.derivative = 1;
%    p.LFnormalize = 0;
%    mdChan = mod(iChan_tmp, 10);
%elseif (iChan_tmp >= 20) && (iChan_tmp < 30)
%    p.LFnormalize = 1;
%    p.derivative = 0;
%    mdChan = iChan_tmp - 20;
%else
%    p.LFnormalize = 0;
%    p.derivative = 0;
%    mdChan = iChan_tmp;
%end

disp(['==== md Channel index: ', num2str(mdChan), ' ===='])
disp('==== Input Parameter: p ====')
disp(p)


%%
%fnameChan = double(p.derivative)*10 + iChan;
%fname0 = ['Chan', num2str(fnameChan)];

fname0 = ['Chan', num2str(iChan)];
disp(fname0)

MDpixelSize_ = MD.pixelSize_;
MDtimeInterval_ = MD.timeInterval_;

% velmap
if mdChan == 0
    
    indPSP = MD.getProcessIndex('ProtrusionSamplingProcess');
    PSP = MD.getProcess(indPSP);
    WSPresult = PSP.loadChannelOutput();
    protmap = WSPresult.avgNormal;
    %size(protmap)
    velmap = protmap*MDpixelSize_/MDtimeInterval_;     %%%%    velocity as nm/sec
    
    disp('size of velocity map')
    size(velmap)
    
    indWP = MD.getProcessIndex('WindowingProcess');
    WP = MD.getProcess(indWP);
    %WP.nBandMax_
    wmax = WP.nSliceMax_;
    tmax = MD.nFrames_;
    
    velmapShift = [nan(wmax, 1), velmap(:, 1:(tmax-1))];
    
    actmap = reshape(velmapShift, wmax, 1, tmax);
    maxLayer = 1;
    
    % actmap
elseif p.WithN == true
    
    % load all channels.mat for allSamplesWithN
    iWinPackInd = MD.getPackageIndex('WindowingPackage');
    if ~isempty(iWinPackInd)
        tmp = MD.packages_{iWinPackInd}.outputDirectory_;
        samplingWithNDirectory_ = fullfile(tmp, 'window_sampling_WithN');
        if ~isdir(samplingWithNDirectory_)
            tmp = MD.outputDirectory_;
            samplingWithNDirectory_ = fullfile(tmp, 'window_sampling_WithN');
        end
    else
        samplingWithNDirectory_ = fullfile(MD.outputDirectory_, 'window_sampling_WithN');
    end
    
    fname = 'all channels.mat';
    inFilePaths = fullfile(samplingWithNDirectory_, fname);
    load(inFilePaths, 'allSamplesWithN');
    
    actmap = allSamplesWithN(mdChan).avg;
    
    disp('size of activity map')
    size(actmap)
    
    indWP = MD.getProcessIndex('WindowingProcess');
    WP = MD.getProcess(indWP);
    %WP.nBandMax_
    wmax = WP.nSliceMax_;
    tmax = MD.nFrames_;
    
else
    
    indWSP = MD.getProcessIndex('WindowSamplingProcess');
    WSP = MD.getProcess(indWSP);
    WSPresult = WSP.loadChannelOutput(mdChan);  % Set channel of  ...
    actmap = WSPresult.avg;                    % actmap > rawActmap, actmap_outl, imActmap
    
    disp('size of activity map')
    size(actmap)
    
    indWP = MD.getProcessIndex('WindowingProcess');
    WP = MD.getProcess(indWP);
    %WP.nBandMax_
    wmax = WP.nSliceMax_;
    tmax = MD.nFrames_;
    
end

%% p.derivative option  (old)

%% Omit windows

if numel(p.omittedWindows) > 0
    actmap(p.omittedWindows, :,:) = NaN;
end


%% Omit time frames

if numel(p.subFrames) > 0
    actmap = actmap(:,:, p.subFrames);
    tmax = size(actmap, 3);             %%% update tmax
end


%% if Folding

if p.Folding == 1
    
    if mod(tmax, 2) == 1
        tmax = tmax + 1;
        [a, b, ~] = size(actmap);
        actmap = cat(3, actmap, nan(a, b));
    end
    foldedMap = nan(wmax, maxLayer, tmax/2);
    
    for w = 1:wmax
        for l = 1:maxLayer
            for t = 1:(tmax/2)
                foldedMap(w, l, t) = mean([actmap(w, l, 2*t-1), actmap(w, l, 2*t)], 'omitnan');
            end
        end
    end
    actmap = foldedMap;
    tmax = size(foldedMap, 3);
    MDtimeInterval_ = MDtimeInterval_ * 2;
    
end



%%  Activity Map Outlier & remove windows
%% omitWin, subFr, Folding -> outlierAdj -> (movAvgSmoothing) -> Diff (smoothing again) -> (impute)
%% 2018/05
%% omitWin, subFr, Folding -> outlierAdj -> (movMedNormalize)  -> (movAvgSmoothing)
%% ((imputed) Common Factor normalization)
%% -> Diff (smoothing again) -> (impute)

%%  outlierAdj
%%  actmap_outl

rawActmap = cell(1, maxLayer);
actmap_outl = cell(1, maxLayer);

for indL = 1:maxLayer
    disp(['==== ', num2str(indL), ' Layer ===='])
    
    rawActmap{indL} = squeeze(actmap(:, indL, :));
    inputmap = rawActmap{indL};
    
    disp('# of NaN in map:')
    disp( sum(isnan(inputmap(:))) )
    
    disp('summary:')
    disp( summary(inputmap(:)) )
    m0 = mean(inputmap(:), 'omitnan'); %disp(m0)
    std0 = std(inputmap(:), 'omitnan'); %disp(std0)
    
    % 5*sigma
    % 2018/05/28. p.outlSigma
    [r, c] = find(abs(inputmap-m0)/std0 > p.outlSigma);
    disp('row, column of outliers:')
    disp([r, c])
    actmap_outl{indL} = inputmap;
    actmap_outl{indL}(abs(inputmap-m0)/std0 > p.outlSigma) = NaN;
    
    disp('upper/lower threshold for outliers:')
    disp(m0+ p.outlSigma*std0)
    disp(m0- p.outlSigma*std0)
    disp('num of outlier')
    disp( sum(sum(abs(inputmap-m0)/std0 > p.outlSigma)) )
    
end



%%  rawActmap, actmap_outl -> normed rawActmap, actmap_outl
% if not p.LFnormalize

LFnormalizingMap = cell(1, maxLayer);

for indL = 1:maxLayer
    LFnormalizingMap{indL} = nan(size(rawActmap{indL}));
end

if (p.LFnormalize == true)
    
    rawActmap_norm = cell(1, maxLayer);
    actmap_outl_norm = cell(1, maxLayer);
    %normalizingFactorMap = cell(1, maxLayer);
    
    for indL = 1:maxLayer
        disp(['==== ', num2str(indL), ' Layer ===='])
        
        inputmap = actmap_outl{indL};
        rawActmap_norm{indL} = inputmap;
        actmap_outl_norm{indL} = inputmap;
        LFnormalizingMap{indL} = nan(size(inputmap));
        
        % movmedian normalization => LowFreqNormalization
        for w = 1:size(inputmap, 1)
            medTS = movmedian(inputmap(w, :), p.movMedFrameSize, 'omitnan');
            LFnormalizingMap{indL}(w, :) = medTS;
        end
        
        % Delta F/F formula (%) 2018/10/26
        rawActmap_norm{indL} = (inputmap - LFnormalizingMap{indL}) ./ LFnormalizingMap{indL} .* 100;
        
        %
        inputmap = rawActmap_norm{indL};
        disp('# of NaN in normed map:')
        disp( sum(isnan(inputmap(:))) )
        
        disp('summary:')
        disp( summary(inputmap(:)) )
        m0 = mean(inputmap(:), 'omitnan'); %disp(m0)
        std0 = std(inputmap(:), 'omitnan'); %disp(std0)
        
        % 5*sigma
        [r, c] = find(abs(inputmap-m0)/std0 > p.outlSigma);
        disp('row, column of outliers:')
        disp([r, c])
        actmap_outl_norm{indL} = inputmap;
        actmap_outl_norm{indL}(abs(inputmap-m0)/std0 > p.outlSigma) = NaN;
        
        disp('upper/lower threshold for outliers:')
        disp(m0+ p.outlSigma*std0)
        disp(m0- p.outlSigma*std0)
        disp('num of outlier')
        disp( sum(sum(abs(inputmap-m0)/std0 > p.outlSigma)) )
        
    end
    
    rawActmap = rawActmap_norm;
    actmap_outl = actmap_outl_norm;
    
end

if p.LFof
    rawActmap = LFnormalizingMap;
    actmap_outl = LFnormalizingMap;
end



%% movAvg smoothing

if (p.movingAvgSmoothing == true)
    actmap_outl2 = actmap_outl;
    
    gaussianFilter = fspecial('gaussian', 3, 0.5);   %% minimal gaussian filtering
    
    for l = 1:maxLayer
        tmp = actmap_outl{l};
        actmap_outl2{l} = nanconv(tmp, gaussianFilter, 'edge', 'nanout');
    end
    actmap_outl = actmap_outl2;
end


%% 10/18/2017 (temp)

if (p.movingAvgSmoothing > 0) & (p.movingAvgSmoothing < 1)
    actmap_outl2 = actmap_outl;
    smParam = p.movingAvgSmoothing;
    
    for indL = 1:maxLayer
        if ~all(isnan(actmap_outl{indL}(:)))
            actmap_outl2{indL} = smoothActivityMap(actmap_outl{indL}, 'SmoothParam', smParam, 'UpSample', 1);
        end
    end
    actmap_outl = actmap_outl2;
end

%% 2019/05/04. EWMA
if (p.movingAvgSmoothing == false) && (p.EWMA < 1)
    actmap_outl2 = actmap_outl;
    lam0 = p.EWMA;
    % EWMA using filter()
    for indL = 1:maxLayer
        tmp = actmap_outl2{indL};
        mvec = mean(tmp, 2, 'omitnan');
        tmp2 = tmp - mvec;
        indmap = isnan(tmp2);
        tmp2(indmap) = 0;
        tmp2Filt = filter(lam0, [1, -(1-lam0)], tmp2, [], 2);
        tmp2Filt(indmap) = NaN;
        tmp2Filt2 = tmp2Filt + mvec;
        actmap_outl{indL} = tmp2Filt2;
    end
end

%%  Ratio map, v3 20190509
if p.Ratio == 1
    % p.baseOfRatio
    [fname3, ~, ~, ~, ~, ~, actmap_outl3, ~, ~, ~, ~] ...
        = mapOutlierImputation(MD, p.baseOfRatio, maxLayer, 'impute', p.impute, ...
        'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
        'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, ...
        'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
        'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
        'EWMA', p.EWMA);
    
    actmap_outlScAdj = cell(1, maxLayer);
    actmap_outlScAdjBase = cell(1, maxLayer);
    actmap_Ratio = cell(1, maxLayer);
    
    tmp = cell2mat(actmap_outl);
    tmp3 = cell2mat(actmap_outl3);
    mNumer = mean(tmp(:), 'omitnan'); sdNumer = std(tmp(:), [], 'omitnan');
    mDenomi = mean(tmp3(:), 'omitnan'); sdDenomi = std(tmp3(:), [], 'omitnan');
    % scaling to Z*10 + 100
    for indL = 1:maxLayer
        actmap_outlScAdj{indL} = (actmap_outl{indL} - mNumer)./sdNumer .* 10 + 100;
        actmap_outlScAdjBase{indL} = (actmap_outl3{indL} - mDenomi)./sdDenomi .* 10 + 100;
        % 9*sigma filtering
        actmap_outlScAdj{indL}(abs(actmap_outlScAdj{indL}-100)./10 > 9) = NaN;
        actmap_outlScAdjBase{indL}(abs(actmap_outlScAdjBase{indL}-100)./10 > 9) = NaN;
        %actmap_Ratio{indL} = actmap_outlScAdj{indL} ./ actmap_outlScAdjBase{indL};
    end
    % use simple regression to make linearly transformed Base signal at
    % global level -> at the deep layer (neg control)
    y = actmap_outlScAdj{maxLayer};
    x = actmap_outlScAdjBase{maxLayer};
    yy = y(:); xx = x(:);
    [bhat, ~, ~] = regress(yy, [ones(size(xx)), xx]);
    % figure output
    figRatio = figure('Visible', 'off');
            %y0 = actmap_outl{maxLayer}; x0 = actmap_outl3{maxLayer};
            %yy0 = y0(:); xx0 = x0(:);
    scatter(xx, yy, 1)
    h = refline([bhat(2), bhat(1)]); h.Color = 'r';
    fnameTmp = ['scaled Chan', num2str(iChan_tmp)];
    layerLab = ['-', num2str(maxLayer), 'L'];
    xlabel(['scaled ', fname3 layerLab]); ylabel([fnameTmp, layerLab]);
    title0 = ['intrcpt: ', num2str(round(bhat(1),1)), ', slope: ', num2str(round(bhat(2), 3))];
    title2 = ['stdy: ', num2str(round(nanstd(yy), 1)), ', stdx: ', num2str(round(nanstd(xx), 1))];
    title({'Ratio metric'; title0; title2})
    set(gca, 'FontSize', 14)
    saveas(figRatio, fullfile(p.figuresDir, 'RatioNormConstants.png'), 'png')
    saveas(figRatio, fullfile(p.figuresDir, 'RatioNormConstants.fig'), 'fig')
    
    % translated X or Base signals
    %yhat = yy - res;
    actmap_outlScAdjHat = cell(1, maxLayer);
    for indL = 1:maxLayer
        actmap_outlScAdjHat{indL} = actmap_outlScAdjBase{indL} .* bhat(2) + bhat(1);
        actmap_Ratio{indL} = actmap_outlScAdj{indL} ./ actmap_outlScAdjHat{indL};
    end
    % figure,imagesc(actmap_outlScAdj{1}),colormap(jet),colorbar
    
    % ratio map output
    actmap_outl = actmap_Ratio;
    % 2019/05/04. EWMA
    %if (p.movingAvgSmoothing == false) && (p.EWMA < 1)
    %    actmap_outl2 = actmap_outl;
    %    lam0 = p.EWMA;
    %    % EWMA using filter()
    %    for indL = 1:maxLayer
    %        tmp = actmap_outl2{indL};
    %        mvec = mean(tmp, 2, 'omitnan');
    %        tmp2 = tmp - mvec;
    %        indmap = isnan(tmp2);
    %        tmp2(indmap) = 0;
    %        tmp2Filt = filter(lam0, [1, -(1-lam0)], tmp2, [], 2);
    %        tmp2Filt(indmap) = NaN;
    %        tmp2Filt2 = tmp2Filt + mvec;
    %        actmap_outl{indL} = tmp2Filt2;
    %    end
    %end
end


%%  2018/05/30  CFnormalize

feigCurve = [];  % empty

if (p.CFnormalize == true)
    
    % movingAvgSmoothing on for CFana if on
    %if (p.movingAvgSmoothing == true)
    %    actmap_outl2 = actmap_outl;
    %    gaussianFilter = fspecial('gaussian', 3, 0.5);   %% minimal gaussian filtering
    %    for l = 1:maxLayer
    %        tmp = actmap_outl{l};
    %        actmap_outl2{l} = nanconv(tmp, gaussianFilter, 'edge', 'nanout');
    %    end
    %    actmap_outl = actmap_outl2;
    %end
    
    % impute always on for CFanal
    imActmap1 = cell(1, maxLayer);
    for indL = 1:maxLayer
        mat = actmap_outl{indL};
        imActmap1{indL} = myknnimpute(mat')';     % Note tha double transeposes are necessary.
    end
    
    % p.CommonFactorNormAddCh
    
%    [fname2, ~, ~, ~, ~, ~, ~, imActmap2] ...
%        = mapOutlierImputation(MD, p.CommonFactorNormAddCh, maxLayer, 'impute', 1, ...
%        'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
%        'Folding', p.Folding, 'subFrames', p.subFrames, ...
%        'movingAvgSmoothing', p.movingAvgSmoothing, 'movMedFrameSize', p.movMedFrameSize, ...
%        'EWMA', p.EWMA);
% 2019/05/11, replaced with ML_mapForCFN() in the pipeline
    S = load(fullfile(p.figuresDir, [p.CommonFactorNormAddCh, '.mat']));
    fname2 = S.fname0;
    imActmap2 = S.imActmap;   
    
    % 2019/05. CFN3 no smoothing in CFN!
    [~, ~, ~, ~, ~, ~, ~, imVelmap] ...
        = mapOutlierImputation(MD, 0, 1, 'impute', 1, ...
        'WithN', 0, 'omittedWindows', p.omittedWindows, ...
        'Folding', p.Folding, 'subFrames', p.subFrames, ...
        'movingAvgSmoothing', p.movingAvgSmoothing, 'EWMA', p.EWMA);
    
    disp(['==== Factor Analysis Starting with the channel, ', fname2])
    
    % FAn function
    [normalizingFactorMap, CFnormMap1, CFnormMap2, feigCurve] = ...
        mapCommonFactorNormalization(imActmap1, imActmap2, imVelmap, wmax, tmax, maxLayer, ...
        'factoranMethod', p.factoranMethod);
    
    %rawActmap = normalizingFactorMap;   % for quick reference
    % 2019/05/09
    rawActmap = CFnormMap1;
    
    actmap_outl = CFnormMap1;
    actmap_outlSctmp = cell(1, maxLayer);
    for indL = 1:maxLayer
        actmap_outlCell = num2cell(actmap_outl{indL}, 2);
        actmap_outlZ = cellfun(@(x) nanZscore(x), actmap_outlCell, 'UniformOutput', false);
        actmap_outlSctmp{indL} = cell2mat(actmap_outlZ);
    end
    
    actmap_outl = actmap_outlSctmp;
end

if p.CFof
    rawActmap = CFnormMap1;     % for quick reference
    actmap_outl = normalizingFactorMap;
end

%% p.VelNorm
VNout = struct();

if (p.VelNorm == true)
        
    % impute always 
    imActmap1 = cell(1, maxLayer);
    for indL = 1:maxLayer
        mat = actmap_outl{indL};
        imActmap1{indL} = myknnimpute(mat')';     % Note tha double transeposes are necessary.
    end

    % 2019/05. 
    [~, ~, ~, ~, ~, ~, ~, imVelmap] ...
        = mapOutlierImputation(MD, 0, 1, 'impute', 1, ...
        'WithN', 0, 'omittedWindows', p.omittedWindows, ...
        'Folding', p.Folding, 'subFrames', p.subFrames, ...
        'movingAvgSmoothing', p.movingAvgSmoothing, 'EWMA', p.EWMA);
    disp(['==== Velocity Normalization for the channel, ', fname0])

    % vel norm with SPAR funtion
    [regAvBICmapcell, regAvBICcurves, BICPvecs, resiMaps] = ...
        mapVelNormalizationWithSPAR_lag0(imActmap1, imVelmap, p.VNtlagMax, wmax, tmax, maxLayer);

    % velmap is regressed out from imActmap1. after/before
    actmap_outl = resiMaps;
    rawActmap = imActmap1;
    
    VNout.regAvBICmapcell = regAvBICmapcell;
    VNout.regAvBICcurves = regAvBICcurves;
    VNout.BICPvecs = BICPvecs;
    
    fBICPvecs = cell(maxLayer, 1);
    for indL = 1:maxLayer
        fBICPvecs{indL} = figure('Visible', 'off');
        tmp1 = squeeze(BICPvecs(indL, :, 1));
        tmp2 = squeeze(BICPvecs(indL, :, 2));
        plot(tmp1); 
        hold on; plot(tmp2)
        title('BIC chosen orders')
        legend('ARorder0', 'ARorder1', 'Location', 'northoutside','Orientation','horizontal')
        saveas(fBICPvecs{indL}, fullfile(p.figuresDir, ['fBICPvecs_', ...
            fname0, '_', num2str(indL),'L.png']), 'png')
    end
end


%% p.derivative option

if (p.derivative == true)
    actmap_outl2 = actmap_outl;
    for l = 1:maxLayer
        tmp = actmap_outl{l};
        difftmp = diff(tmp, [], 2);   % time-wise difference operation
        actmap_outl2{l} = [nan(size(tmp, 1), 1), difftmp];
    end
    actmap_outl = actmap_outl2;
    
    if iChan == 10
        for l = 1:maxLayer
            actmap_outl{l} = actmap_outl{l} ./ MDtimeInterval_;
        end
    end
    
    %%%%%  Smoothing again the differences
    
    if (p.movingAvgSmoothing == true)
        actmap_outl2 = actmap_outl;
        
        gaussianFilter = fspecial('gaussian', 3, 0.5);   %% minimal gaussian filtering
        for l = 1:maxLayer
            tmp = actmap_outl{l};
            actmap_outl2{l} = nanconv(tmp, gaussianFilter, 'edge', 'nanout');
        end
        actmap_outl = actmap_outl2;
    end
    
end



%%  scaling the activity map

actmap_outlSc = cell(1, maxLayer);
for indL = 1:maxLayer
    
    actmap_outlCell = num2cell(actmap_outl{indL}, 2);
    actmap_outlZ = cellfun(@(x) nanZscore(x), actmap_outlCell, 'UniformOutput', false);
    actmap_outlSc{indL} = cell2mat(actmap_outlZ);
    
    % CentMap ...
    %actmap_outlSc{indL} = detrend(actmap_outl{indL}', 'constant')';
end


%%  Imputation (used as an input of computations), later maybe restricted

if p.impute == 1
    
    imActmap = cell(1, maxLayer);
    for indL = 1:maxLayer
        mat = actmap_outl{indL};
        imActmap{indL} = myknnimpute(mat')';     % Note tha double transeposes are necessary.
    end
else
    imActmap = actmap_outl;
end

% if iChan==0, re-define output variables: not necessary

%%  save output
if ~isempty(p.figuresDir) && ~isempty(p.chanName)
    
    if ~isfolder(p.figuresDir); mkdir(p.figuresDir); end
    chanName = p.chanName;
    save(fullfile(p.figuresDir, [chanName, '.mat']), 'fname0', 'MDpixelSize_', 'MDtimeInterval_', 'wmax', ...
        'tmax', 'rawActmap', 'actmap_outl', 'imActmap', 'actmap_outlSc', ...
        'LFnormalizingMap', 'feigCurve', 'iChan', 'maxLayer', 'chanName', 'p', ...
        'VNout')    
end


end

