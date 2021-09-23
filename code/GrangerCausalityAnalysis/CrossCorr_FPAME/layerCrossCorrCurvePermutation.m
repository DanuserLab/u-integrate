function layerCrossCorrCurvePermutation(MD, iChan1, iChan2, chan1Name, chan2Name, layerMax, figuresDir, varargin) 
% layerCrossCorrCurvePermutation() Perform cross correlation analysis between two
% channels at different layers. 
%
% Input:
%       MD          - a movieData object
%       iChan1      - the 1st channel index
%                     iChan = 1x indicates the differenced map (X_t =
%                     X_{t-1}) of channel x.
%                     iChan = 2x indicates movmedian normalized (LFnormalize).
%                     iChan = 3x indicates LowFreq normalizing Map.
%                     iChan = 1xx indicates common factor-based normlized.
%                     iChan = 2xx indicates the CommonFactor normalizing
%                     Map.
%       chan1Name   - a short name for channel1.
%       iChan2      - the 2nd channel index
%       chan2Name   - a short name for channel2.
%       layerMax    - maximum layer to be analyzed      
%       figuresDir  - a directory where plots are saved as png files
%
% Output: png files are saved in the figuresDir. 
%       
% Option:
%       figFlag     - if 'on', matlab figures are ploted. Default is 'off'.
%       lagMax      - the maximum lag to compute. Default is the number of
%                   time frames devided by 4.
%       fullRange   - if true, smoothed xcorr maps are given in [-1, 1]
%                   scale to compare multiple xcorr maps. Default is false.
%       impute      - if true, moderate missing values are imputed by using
%                   knnimpute.m function. Default is false.
%       %parpoolNum  - number of local parallel pool used during permutation. Default is 4.
%       rseed       - input for running rng('default'); rng(rseed). Default
%                   is 'shuffle'. If it is a specific number, the permutation will give
%                   the same result.
%       numPerm     - number of permutation. Default is 1000.
%       WithN       - if true, it uses an alternative windowSampling result
%                   which is obtained by sampleMovieWindowsWithN.m and includes number
%                   of pixels for each windows. Default is false.
%       omittedWindows  
%                   - window index in which activities will be replaced by
%                   NaN. Default is null.
%       subFrames
%                   - specified frames will be only used.        
%
% Updated: 
% J Noh, 2021/02/05. Rename and add comments. 
% J Noh, 2019/04/30. Add 'ratio' chan (4xx) and a 'baseOfRatio' option.
%  
% Jungsik Noh, 2019/03/25

ip = inputParser; 
ip.addParameter('figFlag', 'off');
%ip.addParameter('lagMax', round(tmax/4), @isnumeric);
ip.addParameter('lagMax', 5, @isnumeric);
ip.addParameter('fullRange', false);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('impute', false);
ip.addParameter('WithN', false);
ip.addParameter('subFrames', []);
ip.addParameter('omittedWindows', []);
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('topograph', 'on');
ip.addParameter('Folding', false);
%ip.addParameter('LFnormalize', false);
ip.addParameter('movMedFrameSize', nan);   % movmedian frame size for normalization.
ip.addParameter('outlSigma', 5);
%ip.addParameter('CFnormalize', false);
ip.addParameter('CommonFactorNormAddChVec', {NaN, NaN});
ip.addParameter('factoranMethod', 33);
ip.addParameter('baseOfRatioVec', [NaN, NaN]);  % base chan for ratio comp. default is ch1.


parse(ip, varargin{:})
p = ip.Results;

%figFlag = p.figFlag;

%%  figuresDir setup
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['fieldXcorrCurvePermutation_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')


%%  getting Maps from channels 1, 2

[~, ~,MDtimeInterval_, wmax, tmax, ~, ~, imActmap1] ...
            = mapOutlierImputation(MD, iChan1, layerMax, 'impute', p.impute, ...
            'omittedWindows', p.omittedWindows, 'WithN', p.WithN, 'subFrames', p.subFrames, ...
            'movingAvgSmoothing', p.movingAvgSmoothing, 'Folding', p.Folding, ...
            'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
            'CommonFactorNormAddCh', p.CommonFactorNormAddChVec{1}, 'factoranMethod', p.factoranMethod, ...
            'baseOfRatio', p.baseOfRatioVec(1));    

[~, ~, ~, ~, ~, ~, ~, imActmap2] ...
            = mapOutlierImputation(MD, iChan2, layerMax, 'impute', p.impute, ...
            'omittedWindows', p.omittedWindows, 'WithN', p.WithN, 'subFrames', p.subFrames, ...
            'movingAvgSmoothing', p.movingAvgSmoothing, 'Folding', p.Folding, ...
            'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
            'CommonFactorNormAddCh', p.CommonFactorNormAddChVec{2}, 'factoranMethod', p.factoranMethod, ...
            'baseOfRatio', p.baseOfRatioVec(2));     

% 
if p.lagMax == 5    % only when not specified
    p.lagMax = round(tmax/4);
end        

%% variable set up

%imActmap1_layer = cell(1, layerMax);
%imActmap2_layer = cell(1, layerMax);

%for indL = 1:layerMax
%    imActmap1_layer{indL} = reshape(imActmap1{indL}, wmax, 1, tmax);
%end
%imActmap1_3dim = cat(2, imActmap1_layer{1:end});

%for indL = 1:layerMax
%    imActmap2_layer{indL} = reshape(imActmap2{indL}, wmax, 1, tmax);
%end
%imActmap2_3dim = cat(2, imActmap2_layer{1:end});
%%  to handle vel channel reads

if layerMax > 1
    if numel(imActmap1) == 1
        for indL = 2:layerMax; imActmap1{indL} = imActmap1{1}; end
    end
    if numel(imActmap2) == 1
        for indL = 2:layerMax; imActmap2{indL} = imActmap2{1}; end
    end
end


%%  input prepare and call xcorrCurvePermutationTest

ch1Actmap = imActmap1;
ch2Actmap = imActmap2;
ch1ActmapName = chan1Name;
ch2ActmapName = chan2Name;
%lagMax
fsaveName0 = ['Ch', num2str(iChan1), 'Ch', num2str(iChan2)];

%xcorrMat = xcorrCurvePermutationTest(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
xcorrMatArr = layerCrossCorrCurves(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
              fsaveName0, MDtimeInterval_, figuresDir, ...
              'figFlag', p.figFlag, 'fullRange', p.fullRange, 'lagMax', p.lagMax, ...
              'numPerm', p.numPerm, 'parpoolNum', p.parpoolNum, 'rseed', p.rseed);
  


%%  Topographs of xcorr
if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);

nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
topoMap = nan(wmax, nBandMax_);

for indL = 1:layerMax
    tmp = xcorrMatArr{indL, indL};
    tmp2 = tmp(:, p.lagMax+1);                % 7 = lagMax+1+h. lag=gef_t+..., stdN_t
    topoMap(:, indL) = tmp2;
end


title0 = ['xcorr(', chan1Name, '_{t}, ', chan2Name, '_t)'];   % lag = {t+ ...} 
topoFig_xcorrChan1Chan2 = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);
 

%%
saveas(topoFig_xcorrChan1Chan2, fullfile(figuresDir, ['/topoFig_xcorr_', fsaveName0, '.png']), 'png')  

end

%%
disp('====End of layerCrossCorrCurvePermutation ====')


end




