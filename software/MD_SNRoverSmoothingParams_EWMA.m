function MD_SNRoverSmoothingParams_EWMA(MD, iChan, maxLayer, chanName, arOrder, figuresDir, ...
    varargin)
% MD_SNRoverSmoothingParams_EWMA() Compute Singal-to-Noise-Ratios (SNRs) 
% of a given time series map based on Auto-Regressive modeling. 
% Multiple SNRs are for different EWMA smoothing parameters from 1 (no smoothing) 
% to 0.1 (heavy sm) by default. 
%
% Updates:
% J Noh, 2021/01/27. Rename and add comment.
% J Noh, 2019/04/18. Remove Gaussian filter 3by3 smoothing case. Add trend
% line.
% Jungsik Noh, 2019/04/16.
%
% Copyright (C) 2023, Danuser Lab - UTSouthwestern 
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

ip = inputParser;
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
ip.addParameter('adf', 0);
ip.addParameter('figFlag', 'off');
ip.addParameter('parpoolNum', 4);
ip.addParameter('rseed', 'shuffle');
ip.addParameter('numPerm', 1000);
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('smParam', 0.8);
ip.addParameter('LFnormalize', false);
ip.addParameter('movMedFrameSize', nan);   % movmedian frame size for normalization.
ip.addParameter('outlSigma', 5);
ip.addParameter('CFnormalize', false);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('factoranMethod', 33);
ip.addParameter('smParVec', 1:-0.05:0.1);   % 31 pts
ip.addParameter('baseOfRatio', nan);  % base chan for ratio comp. default is ch1.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 

ip.parse(varargin{:});
p = ip.Results;

%figFlag = p.figFlag;
set(groot,'defaultLegendAutoUpdate','off')

%%  figuresDir setup

if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['MD_SNRoverSmoothingParams_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')

%%  getting Maps from channels

disp(chanName)

[fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap, ~, lowFreqmap, feigCurve] ...
    = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
    'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
    'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing', true, ...
    'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio,  ...
            'figuresDir', figuresDir, 'chanName', chanName); 

disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])

%%  Remove columns at the beginning or end with all NaN's (global nanK, nanK2)

nanK = 0;

for indL = 1:maxLayer
    map = imActmap{indL};
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK = max(nanK, map_ik);
end

disp(['nanK: ', num2str(nanK)])

%
nanK2 = 0;

for indL = 1:maxLayer
    map = flip(imActmap{indL}, 2);
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK2 = max(nanK2, map_ik);
end

disp(['nanK2: ', num2str(nanK2)])

%%  getting maps with all kinds of smoothing params

numSms = numel(p.smParVec) ;
smImActMaps = cell(maxLayer, numSms);

[~, ~, ~, ~, ~, ~, ~, imActmap_sm1, ~, ~, ~] ...
    = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
    'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
    'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing', false, ...
    'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio, ...
            'figuresDir', figuresDir, 'chanName', chanName); 
        
for indL = 1:maxLayer
    smImActMaps{indL, 1} = imActmap_sm1{indL}(:, nanK+1:tmax-nanK2);
end

%% 2:numSms with smParam values

for k = 2:numSms
    
    [~, ~, ~, ~, ~, ~, ~, imActmap_tmp, ~, ~, ~] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
        'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
        'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing',false,'EWMA', p.smParVec(k), ...
        'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
        'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatio, ...
            'figuresDir', figuresDir, 'chanName', chanName); 
    
    for indL = 1:maxLayer
        smImActMaps{indL, k} = imActmap_tmp{indL}(:, nanK+1:tmax-nanK2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  SNRoverSmoothParam
%arOrder = 32

SNRmat = cell(maxLayer, 1);

for indL = 1:maxLayer
    for sm = 1:numSms
        smImActMaps_tmp = smImActMaps{indL, sm};
        SNRmat_tmp = nan(wmax, 1);
        
        parfor w = 1:wmax
            
            y = smImActMaps_tmp(w, :)';
            N = numel(y);
            if any(isnan(y))
                snr0 = nan;
            else
                yLag = myLagmatrix(y, arOrder);
                X = [ones(N,1), yLag];
                bhat = X\y;
                sse = (y-X*bhat)'*(y-X*bhat);
                errVar = sse/(N-size(X,2));
                yhat = X*bhat;
                % Signal-to-noise ratio = std(X*beta) / RMSE
                snr0 = std(yhat)/sqrt(errVar);
                %Rsq0(1) = 1 - sse/( (y-mean(y))'*(y-mean(y)) );
                %avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
                %avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
                %avAIC(1) = -2*avlogLhat + 2/N*(size(X,2) + 1);
            end
            
            SNRmat_tmp(w) = snr0;
        end
        SNRmat{indL}(:, sm) = SNRmat_tmp;
        %disp(sm)
    end
end

%%  BoxPlot per layer

fbpsnr = cell(maxLayer, 1);
% penalization on roughness
xticklab = [1 - p.smParVec];
medSNR = cell(maxLayer, 1);
tmp = cell2mat(SNRmat);
ymax = log10(max(tmp(:)));
ymin = log10(min(tmp(:)));

for indL = 1:maxLayer
    medSNR{indL} = reshape(median(SNRmat{indL}, 1, 'omitnan'), [], 1);
    
    fbpsnr{indL} = figure('Visible', p.figFlag);
    boxplot(log10(SNRmat{indL}))
    ylim([ymin, ymax])
    
    hold on
    plot(log10(medSNR{indL}), 'r')
    refline([0,0])
    % snr = 2 or 1/2
    h1 = refline([0, log10(2)]); h1.LineStyle = '--';
    h2 = refline([0, log10(1/2)]); h2.LineStyle = '--';
    
    set(gca, 'XTickLabelRotation', 45)
    set(gca, 'XTickLabel', xticklab)
    set(gca, 'FontSize', 12)
    
    title([chanName, '-', num2str(indL), 'L', ' ARorder: ', num2str(arOrder)])
    xlabel('Smoothness (1-lambda)')
    ylabel('log10(SNR) (std(Sig)/std(err))')
    %
    lsnr = log10(medSNR{indL}); lsnr1 = lsnr(end-3:end);
    m0 = (lsnr1(end) - lsnr1(1))/3;
    tmp = lsnr(end) - m0 * (0:numSms-1);
    linTrend = flip(tmp);
    plot(linTrend, 'k--');
    
    pause(1)
    
end

%%
for indL = 1:maxLayer
    saveas(fbpsnr{indL}, fullfile(figuresDir, [fname0, '_SNRoverSmoothParam_', num2str(indL), 'L.png']), 'png')
end


%% save
smParVec = p.smParVec;
save(fullfile(figuresDir, [fname0, '_medSNR.mat']), 'medSNR', 'smParVec', 'SNRmat')

%% SNR over windows
% SNRmat

%fsnrWind = cell(maxLayer, 1);
fsnrWind = figure('Visible', p.figFlag);

for indL = 1:maxLayer
    snrVec = SNRmat{indL}(:, 1);
    plot(1:wmax, snrVec)
    hold on
end

title([chanName, ', ARorder: ', num2str(arOrder)])
a = 1:maxLayer; b = num2str(a(:)); c = cellstr(b); d = strcat(c, {'L'});
legend(d, 'Location','northoutside','Orientation','horizontal');
ax = gca;
ax.FontSize = 15;

xlabel('Window'); ylabel('SNR')
refline([0,1])
h1 = refline([0, (2)]); h1.LineStyle = '--';
h2 = refline([0, (1/2)]); h2.LineStyle = '--';

%
saveas(fsnrWind, fullfile(figuresDir, [fname0, '_fsnrWind.png']), 'png')

%%
disp('====End of MD_SNRoverSmoothingParams_EWMA ====')

end


