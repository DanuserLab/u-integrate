function mapXcorrCurvePermutation_Vel(MD, iChan1, chan1Name, layerMax, figuresDir, varargin) 
% mapXcorrCurvePermutation_Vel Perform cross correlation analysis between a
% channel and the edge velocity. It plots cross correlation maps, their mean curves
% at each layer together with confidence bounds based on permutation, and a
% topograph of the cross correlations at lag 0. The cross correlations at
% lag h are Corr(chan1_{t+h}, Vel_t).
% It computes cross correlations in a fashion that can handle many NaN's 
% by utilizing nanXcorrMaps.m function.
%
% Usage:
%       mapXcorrCurvePermutation_Vel(MD, 1, 'Actin', 4, figuresDir, 'impute', 0)
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
%       chan1Name   - a short name for the channel. eg. 'Actin'
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
%                   knnimpute.m function. Default is true.
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
%       topograph   - if 'off' topographs are not plotted. Default is 'on'.
%
% Updated: 
% J. Noh, 2019/05/31. Add a 'VNtlagMax' option for Vel-Normalization.
% J Noh, 2019/05/04. 'EWMA' specifies lambda value for exponentially
% weighted moving average smoothing. Default is 1 (raw data).
% J Noh, 2019/04/30. Add 'ratio' chan (4xx) and a 'baseOfRatio' option.
% J Noh, 2019/02/28. Add an option, 'factoranMethod'.
% J Noh, 2018/11/14. Remove parpoolNum, parpool().
% J Noh, 2018/10/28. Move 'LFnormalize', 'CFnormalize' options inside of
% mapOutlierImputation().
% J Noh, 2018/10/26. Add a 'LFnormalize' option. iChan includes LFof,
% CFof.
% J Noh, 2018/05/30. Add Common Factor (CF) analysis-based normalization.
% J Noh, 2018/05/28. Add 'outlSigma' option.
% J Noh, 2018/05/27. Introduce 'normalization' and p.movMedFrameSize.
% J Noh, 2018/01/14. 'lagMax' is adjusted for the case when 'subFrames' is
% specified.
% J Noh, 2017/10/11, raw activities can be smoothed. New option is
% 'movingAvgSmoothing'.
% J Noh, 2017/08/26. To deal with differenced channels. 
%                     iChan = 1x indicates the differenced map (X_t =
%                     X_{t-1}) of channel x.
% Jungsik Noh, 2017/05/23  
% Jungsik Noh, 2016/10/22
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
%
% This file is part of GrangerCausalityAnalysisPackage.
% 
% GrangerCausalityAnalysisPackage is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GrangerCausalityAnalysisPackage is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GrangerCausalityAnalysisPackage.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

%tmax = MD.nFrames_;

ip = inputParser; 
ip.addParameter('figFlag', 'off');
%ip.addParameter('lagMax', round(tmax/4), @isnumeric);
ip.addParameter('lagMax', 5, @isnumeric);
ip.addParameter('fullRange', false);
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
ip.addParameter('omittedWindows', []);
ip.addParameter('h0', []);
ip.addParameter('mvFrSize', 0);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', false);
%ip.addParameter('LFnormalize', false);
ip.addParameter('movMedFrameSize', nan);   % movmedian frame size for normalization.
ip.addParameter('outlSigma', 5);
%ip.addParameter('CFnormalize', false);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('factoranMethod', 33);
ip.addParameter('baseOfRatio', 1);  % base chan for ratio comp. default is ch1.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 
ip.addParameter('VNtlagMax', 10);   % maximum lag (frames) for VelNorm with SPAR


parse(ip, varargin{:})
p = ip.Results;

%figFlag = p.figFlag;

%%  figuresDir setup
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['mapXcorrCurvePermutation_Vel_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')


%%  getting Maps from channel & vel (ch0)

[~, ~,MDtimeInterval_, wmax, tmax, ~, ~, imActmap1] ...
            = mapOutlierImputation(MD, iChan1, layerMax, 'impute', p.impute, 'WithN', p.WithN, ...
                'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
                'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, ...
                'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
                'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
            'baseOfRatio', p.baseOfRatio, 'EWMA', p.EWMA, ...
            'figuresDir', figuresDir, 'chanName', chan1Name, 'VNtlagMax', p.VNtlagMax); 

[~, ~, ~, ~, ~, ~, ~, imVelmap] ...
            = mapOutlierImputation(MD, 0, 1, 'impute', p.impute, 'omittedWindows', ...
            p.omittedWindows, 'Folding', p.Folding, 'subFrames', p.subFrames, ...
            'movingAvgSmoothing', p.movingAvgSmoothing, 'outlSigma', p.outlSigma, 'EWMA', p.EWMA);    


% 
if p.lagMax == 5    % only when not specified
    p.lagMax = round(tmax/4);
end

        
        
%%  adjust actmap according to velmap. Match layers

for indL = 1:layerMax
    imActmap1{indL}(:, 1) = [];             % Start from frame=2
end

imActmap2 = cell(1, layerMax);
for indL = 1:layerMax
    imActmap2{indL} = imVelmap{1}(:, 2:tmax);
end
      

%%  input prepare and call xcorrCurvePermutationTest

ch1Actmap = imActmap1;
ch2Actmap = imActmap2;
ch1ActmapName = chan1Name;
ch2ActmapName = 'Vel';
%lagMax
fsaveName0 = ['Ch', num2str(iChan1), 'Ch', num2str(0)];

xcorrMat = xcorrCurvePermutationTest(ch1Actmap, ch2Actmap, ch1ActmapName, ch2ActmapName, ...
              fsaveName0, MDtimeInterval_, figuresDir, ...
              'figFlag', p.figFlag, 'fullRange', p.fullRange, 'lagMax', p.lagMax, ...
              'numPerm', p.numPerm, 'parpoolNum', p.parpoolNum, 'rseed', p.rseed);
  


%%  Topographs of xcorr
if strcmp(p.topograph, 'on')

    iWinProc = MD.getProcessIndex('WindowingProcess',1,0);

    nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
    topoMap = nan(wmax, nBandMax_);

    for indL = 1:layerMax
        tmp = xcorrMat{indL};
        tmp2 = tmp(:, p.lagMax+1);                % 7 = lagMax+1+h. lag=gef_t+..., stdN_t
        topoMap(:, indL) = tmp2;
    end


    title0 = ['xcorr(', ch1ActmapName, '_{t}, ', ch2ActmapName, '_t)'];   % lag = {t+ ...} 
    topoFig_xcorrChan1Chan2 = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);


    %%
    saveas(topoFig_xcorrChan1Chan2, fullfile(figuresDir, ['/topoFig_xcorr_', fsaveName0, '.png']), 'png')  

end

%%
disp('====End of mapXcorrCurvePermutation_Vel====')


end




