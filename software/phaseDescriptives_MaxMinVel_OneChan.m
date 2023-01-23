function phaseDescriptives_MaxMinVel_OneChan(MD, iChan, maxLayer, chanName, ...
            chanTitle, smParamTh, samplingBw, figuresDir, varargin)
% phaseDescriptives_MaxMinVel_OneChan Draw descriptive plots of profiled
% kinetics around the timing of maximum or minimum velocity.
%
% Options:  subFrames
%                   - specified frames will be only used.        
%
% Updated: 
% J Noh, 2021/07/21. Due to Java memory error in GUI, add figure close
% commands.
% J Noh, 2019/04/30. Add 'ratio' chan (4xx) and a 'baseOfRatio' option.
% J Noh, 2019/04/18. Add 'movingAvgSmoothing' option in phaseMaskinginternal() 
% too  to take smoothed activity map input. 
% J Noh, 2019/04/03. save criptsMask2, add minimumRL (in ppline too), and
% heatmaps for ch1Zsamples. 
% J Noh, 2019/02/28. Add an option, 'factoranMethod'.
% J Noh, 2018/10/29. Add options of 'CommonFactorNormAddCh', 'outlSigma'. 
% J Noh, 2018/05/27. Introduce 'normalization' and p.movMedFrameSize.
% J Noh, 2018/04/03. Fix all nan cases.
% Jungsik Noh, 2018/02/05. Instead of truncating samplingBw time frames at
% the beginning and end, NaN's with samplingBw length are added front and
% reat, so that we can sample all the motion events.
% Jungsik Noh, 2018/01/30. Figure output are improved. A few figures have .fig format.
% Jungsik Noh, 2017/05/23
% Jungsik Noh, 2017/04/05
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
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
%ip.addParameter('parpoolNum', 4);
%ip.addParameter('rseed', 'shuffle');
%ip.addParameter('numPerm', 1000);

ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('movMedFrameSize', nan);   % movmedian frame size for normalization.
ip.addParameter('outlSigma', 5);
ip.addParameter('CommonFactorNormAddCh', NaN);
ip.addParameter('factoranMethod', 33);
ip.addParameter('minimumRunLength', 5);  
ip.addParameter('baseOfRatio', 1);  % base chan for ratio comp. default is ch1.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 


ip.parse(varargin{:});
p = ip.Results;

figFlag = p.figFlag;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['phaseDescriptives_MaxMinVel_OneChan_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')


%%  getting Maps from channels

%
disp(chanName)
disp(chanTitle)


[fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap, actmap_outlSc] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, 'WithN', p.WithN, ...
                'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
                'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, ...
                'movMedFrameSize', p.movMedFrameSize, 'outlSigma', p.outlSigma, ...
            'CommonFactorNormAddCh', p.CommonFactorNormAddCh, 'factoranMethod', p.factoranMethod, ...
            'baseOfRatio', p.baseOfRatio, 'EWMA', p.EWMA);  
            
disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])

%
[~, ~, ~, ~, ~, ~, ~, imVelmap, velmap_outlSc] ...
            = mapOutlierImputation(MD, 0, 1, 'impute', p.impute, ...
            'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
            'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, 'EWMA', p.EWMA); 

%imvelocitymap = imVelmap{1}(:, 2:tmax);  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 


%%  phaseMaskingInternal

[ProtMask, RetMask, smVelmap] = phaseMaskingInternal(MD, smParamTh, 'minimumRunLength', p.minimumRunLength, ...
    'omittedWindows', p.omittedWindows, 'Folding', p.Folding, 'subFrames', p.subFrames, ...
    'movingAvgSmoothing', p.movingAvgSmoothing, 'EWMA', p.EWMA); 

%%  Prot Onset/Max Vel analysis

%onsetOfMask = [nan(wmax, 1), diff(Mask, [], 2)];

% ProtMask
tmpMask = ProtMask; tmpMask(:, 1) = NaN;
onsetOfMask = [nan(wmax, 1), diff(tmpMask, [], 2)];

%
onset2 = onsetOfMask;
c = zeros(wmax, 1);
d = zeros(wmax, 1);
for w=1:wmax
    xx = onsetOfMask(w, :);
    if isempty( find( (xx == 1), 1) )
        onset2(w, :) = 0;
    else
        c(w) = find( (xx == 1), 1);
        onset2(w, 1:(c(w)-1)) = 0;
    end
    
    if isempty( find( (xx == -1), 1, 'last') )
        onset2(w, :) = 0;
    else
        d(w) = find( (xx == -1), 1, 'last');
        onset2(w, (d(w)+1):tmax) = 0;
    end
end

%
mask2 = cumsum(onset2, 2);

%figure, imagesc(mask2), axis xy, colormap(jet)

timeSegments = cell(wmax, 1);
criptsMask = zeros(size(onsetOfMask));
for w = 1:wmax
    xx = mask2(w, :);
    xxrle = rle(xx);
    val = xxrle(1:2:end);
    vallength = xxrle(2:2:end);
    
    numSegments = sum(val == 1);
    if numSegments == 0
        timeSegments{w}{1} = NaN;
    else
        for k = 1:numSegments
            startFr_k = sum(vallength(1:(2*k-1))) + 1;
            endFr_k = sum(vallength(1:(2*k)));
            timeSegments{w}{k} = [startFr_k, endFr_k];
            
            velSeg = smVelmap(w, startFr_k:endFr_k);
            [~, indFr] = max(velSeg);               % local max in smoothed velSeg
            if (indFr > 1) && (indFr < numel(velSeg))
                fr0 = startFr_k + indFr - 1;
                criptsMask(w, fr0) = 1;
            end
        end
    end
end


%% Mask plot

criptsMask2 = onset2 + criptsMask.*0.5;
%figure, imagesc(tmp), axis xy, colormap(jet)

inputmap = criptsMask2;
fprotMaxMask= figure('Visible', p.figFlag);
figtmp = imagesc(inputmap);
title(['Protrusion onsets / Time of max Vel -', chanName])
axis xy;xlabel('Time (s)');ylabel('Window')
%
cmap0 = jet(41);
colorbar;colormap(cmap0([1, 11, 21, 31, 41], :))

figtmp.AlphaData = 1-isnan(inputmap);

ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%% save criptsMask2
save(fullfile(figuresDir, ['protMaxVel_criptsMask2', chanName, '.mat']), 'criptsMask2')
saveas(fprotMaxMask, fullfile(figuresDir, ['protMaxVelMask_', chanName, '.png']), 'png')
saveas(fprotMaxMask, fullfile(figuresDir, ['protMaxVelMask_', chanName, '.fig']), 'fig')

%% close figures
pause(0.5) 
close(fprotMaxMask) 

%%  MaxMin Vel analysis: (2) sampling Zscores around the onset
% moving band width
%samplingBw = 8;

ch0Zsamples = cell(maxLayer, 1);
ch1Zsamples = cell(maxLayer, 1);

[r, c] = find(criptsMask == 1);
%ind1 = ((c > samplingBw) & (c <= tmax - samplingBw));
%r = r(ind1);
%c = c(ind1);


for indL = 1:maxLayer
    %
    velmap_outlSc_aug = [nan(wmax, samplingBw), velmap_outlSc{1}, ...
        nan(wmax, samplingBw)];
    actmap_outlSc_aug = [nan(wmax, samplingBw), actmap_outlSc{indL}, ...
        nan(wmax, samplingBw)];
    
    ch0Zsamples{indL} = cell(numel(r), 1);
    ch1Zsamples{indL} = cell(numel(r), 1);
    for k = 1:numel(r)
        %timeInt0 = c(k)-samplingBw;
        %timeInt1 = c(k)+samplingBw;
        timeInt0 = c(k);
        timeInt1 = c(k) + 2 * samplingBw;
        ch0Zsamples{indL}{k} = velmap_outlSc_aug(r(k), timeInt0:timeInt1);
        ch1Zsamples{indL}{k} = actmap_outlSc_aug(r(k), timeInt0:timeInt1);
    end
end


%% Onset analysis: plot sampled Zscores 

timeAxis = (-samplingBw:1:samplingBw)*MDtimeInterval_;


for indL = 1:maxLayer
    
    mat1 = cell2mat(ch0Zsamples{indL});
    mat2 = cell2mat(ch1Zsamples{indL});

    meanTS_max_fig{indL} = figure('Visible', p.figFlag); 
    
    if ~all(isnan([mat1; mat2]))
        mat1tmp = mat1(~any(isnan(mat1), 2), :);
        mat2tmp = mat2(~any(isnan(mat2), 2), :);
        lrstd1 = sqrt(diag(long_run_variance(mat1tmp)));
        lrstd2 = sqrt(diag(long_run_variance(mat2tmp)));

        % all nan case
        if isempty(lrstd1) | isempty(lrstd2)
            lrstd1 = zeros(1, size(mat1, 2)); lrstd2 = zeros(1, size(mat2, 2));
        end

        s1 = shadedErrorBarV2(timeAxis, mean(mat1, 1, 'omitnan'), 2*lrstd1/sqrt(size(mat1, 1)), ...
            'lineprops', '-r');
        hold on
        s2 = shadedErrorBarV2(timeAxis, mean(mat2, 1, 'omitnan'), 2*lrstd2/sqrt(size(mat2, 1)), ...
            'lineprops', '-b');

        title1 = ['Average of locally sampled standardized TS'];
        title0 = [num2str(indL), 'L, bw (lag) around maximum velocity time: ', num2str(samplingBw)];
        title({title1; title0})

        refline([0, 0])
        hold on; h = line([0 0], ylim);
        xlabel('Time (s)');ylabel('Standardized TS')
        legend([s1.mainLine, s2.mainLine], {'Vel', chanName})
    end

end


%% saveas
for indL = 1:maxLayer
    saveas(meanTS_max_fig{indL}, fullfile(figuresDir, [chanName, '-maxVel_meanTS_',  num2str(indL), 'L.png']), 'png')
end

save(fullfile(figuresDir, [chanName, '-maxVel-ch0Zsamples.mat']), 'ch0Zsamples');    
save(fullfile(figuresDir, [chanName, '-maxVel-ch1Zsamples.mat']), 'ch1Zsamples');

%% close figures
pause(0.5)
for indL = 1:maxLayer
    close(meanTS_max_fig{indL})
end

%% plot ch1Zsamples heatmap, 2019/04/03
%
smParam = 0.9;

lagMax = samplingBw;
lagGrid = floor(samplingBw/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeInterval_, 2);

fchanZ = cell(1, maxLayer);

for indL = 1:maxLayer

    mat2 = cell2mat(ch1Zsamples{indL});
    inputmap = mat2;
    %alphadata0 = ones(size(inputmap));
    %alphadata0(~Mask) = 0.2;
    %alphadata0(isnan(inputmap)) = 0;

    % nan case
    if all(isnan(inputmap))
        smmap = nan(size(inputmap));
    else
        smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
    end

    fchanZ{indL} = figure('Visible', p.figFlag);
    figtmp = imagesc(smmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;

    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    set(gca, 'XGrid', 'on')

    ax.FontSize = 15;

end

%%
for indL = 1:maxLayer
    saveas(fchanZ{indL}, fullfile(figuresDir, [chanName, '_maxVel_ch1Zsamples', num2str(indL), 'L.png']), 'png')
    saveas(fchanZ{indL}, fullfile(figuresDir, [chanName, '_maxVel_ch1Zsamples', num2str(indL), 'L.fig']), 'fig')
end

%% close figures
pause(0.5)
for indL = 1:maxLayer
    close(fchanZ{indL})
end

%%  Ret Onset/Min Vel analysis

%onsetOfMask = [nan(wmax, 1), diff(Mask, [], 2)];

% RetMask
tmpMask = RetMask; tmpMask(:, 1) = NaN;
onsetOfMask = [nan(wmax, 1), diff(tmpMask, [], 2)];

%
onset2 = onsetOfMask;
c = zeros(wmax, 1);
d = zeros(wmax, 1);
for w=1:wmax
    xx = onsetOfMask(w, :);
    if isempty( find( (xx == 1), 1) )
        onset2(w, :) = 0;
    else
        c(w) = find( (xx == 1), 1);
        onset2(w, 1:(c(w)-1)) = 0;
    end
    
    if isempty( find( (xx == -1), 1, 'last') )
        onset2(w, :) = 0;
    else
        d(w) = find( (xx == -1), 1, 'last');
        onset2(w, (d(w)+1):tmax) = 0;
    end
end

%
mask2 = cumsum(onset2, 2);

%figure, imagesc(mask2), axis xy, colormap(jet)

timeSegments = cell(wmax, 1);
criptsMask = zeros(size(onsetOfMask));
for w = 1:wmax
    xx = mask2(w, :);
    xxrle = rle(xx);
    val = xxrle(1:2:end);
    vallength = xxrle(2:2:end);
    
    numSegments = sum(val == 1);
    if numSegments == 0
        timeSegments{w}{1} = NaN;
    else
        for k = 1:numSegments
            startFr_k = sum(vallength(1:(2*k-1))) + 1;
            endFr_k = sum(vallength(1:(2*k)));
            timeSegments{w}{k} = [startFr_k, endFr_k];
            
            velSeg = smVelmap(w, startFr_k:endFr_k);
            [~, indFr] = min(velSeg);                   % minimum vel
            if (indFr > 1) && (indFr < numel(velSeg))
                fr0 = startFr_k + indFr - 1;
                criptsMask(w, fr0) = 1;
            end
        end
    end
end


%% Mask plot

criptsMask2 = onset2 + criptsMask.* -0.5;
%figure, imagesc(tmp), axis xy, colormap(jet)

inputmap = criptsMask2;
fretMinMask= figure('Visible', p.figFlag);
figtmp = imagesc(inputmap);
title(['Retraction onsets / Time of min Vel -', chanName])
axis xy;xlabel('Time (s)');ylabel('Window')
%
cmap0 = jet(41);
colorbar;colormap(cmap0([1, 11, 21, 31, 41], :))

figtmp.AlphaData = 1-isnan(inputmap);

ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%% save criptsMask2
save(fullfile(figuresDir, ['retMinVel_criptsMask2', chanName, '.mat']), 'criptsMask2')
saveas(fretMinMask, fullfile(figuresDir, ['retMinVelMask_', chanName, '.png']), 'png')
saveas(fretMinMask, fullfile(figuresDir, ['retMinVelMask_', chanName, '.fig']), 'fig')

%% close figures
pause(0.5) 
close(fretMinMask) 

%%  MaxMin Vel analysis: (2) sampling Zscores around the onset
% moving band width
%samplingBw = 8;

ch0Zsamples = cell(maxLayer, 1);
ch1Zsamples = cell(maxLayer, 1);

[r, c] = find(criptsMask == 1);
%ind1 = ((c > samplingBw) & (c <= tmax - samplingBw));
%r = r(ind1);
%c = c(ind1);

for indL = 1:maxLayer
    %
    velmap_outlSc_aug = [nan(wmax, samplingBw), velmap_outlSc{1}, ...
        nan(wmax, samplingBw)];
    actmap_outlSc_aug = [nan(wmax, samplingBw), actmap_outlSc{indL}, ...
        nan(wmax, samplingBw)];
    
    ch0Zsamples{indL} = cell(numel(r), 1);
    ch1Zsamples{indL} = cell(numel(r), 1);
    for k = 1:numel(r)
        %timeInt0 = c(k)-samplingBw;
        %timeInt1 = c(k)+samplingBw;
        timeInt0 = c(k);
        timeInt1 = c(k) + 2 * samplingBw;
        ch0Zsamples{indL}{k} = velmap_outlSc_aug(r(k), timeInt0:timeInt1);
        ch1Zsamples{indL}{k} = actmap_outlSc_aug(r(k), timeInt0:timeInt1);
    end
end



%% Onset analysis: plot sampled Zscores 

timeAxis = (-samplingBw:1:samplingBw)*MDtimeInterval_;


for indL = 1:maxLayer
    
    mat1 = cell2mat(ch0Zsamples{indL});
    mat2 = cell2mat(ch1Zsamples{indL});
    
    meanTS_min_fig{indL} = figure('Visible', p.figFlag); 
    
    if ~all(isnan([mat1; mat2]))
        mat1tmp = mat1(~any(isnan(mat1), 2), :);
        mat2tmp = mat2(~any(isnan(mat2), 2), :);
        lrstd1 = sqrt(diag(long_run_variance(mat1tmp)));
        lrstd2 = sqrt(diag(long_run_variance(mat2tmp)));

        %meanTS_min_fig{indL} = figure('Visible', p.figFlag); 
        s1 = shadedErrorBarV2(timeAxis, mean(mat1, 1, 'omitnan'), 2*lrstd1/sqrt(size(mat1, 1)), ...
            'lineprops', '-r');
        hold on
        s2 = shadedErrorBarV2(timeAxis, mean(mat2, 1, 'omitnan'), 2*lrstd2/sqrt(size(mat2, 1)), ...
            'lineprops', '-b');

        title1 = ['Average of locally sampled standardized TS'];
        title0 = [num2str(indL), 'L, bw (lag) around minimum velocity time: ', num2str(samplingBw)];
        title({title1; title0})

        refline([0, 0])
        hold on; h = line([0 0], ylim);
        xlabel('Time (s)');ylabel('Standardized TS')
        legend([s1.mainLine, s2.mainLine], {'Vel', chanName})
    end

end


%% saveas
for indL = 1:maxLayer
    saveas(meanTS_min_fig{indL}, fullfile(figuresDir, [chanName, '-minVel_meanTS_',  num2str(indL), 'L.png']), 'png')
end

save(fullfile(figuresDir, [chanName, '-minVel-ch0Zsamples.mat']), 'ch0Zsamples');    
save(fullfile(figuresDir, [chanName, '-minVel-ch1Zsamples.mat']), 'ch1Zsamples');

%% close figures
pause(0.5)
for indL = 1:maxLayer
    close(meanTS_min_fig{indL})
end

%% plot ch1Zsamples heatmap, 2019/04/03
%
smParam = 0.9;

lagMax = samplingBw;
lagGrid = floor(samplingBw/2);
xcmapXtick = [1, 1+lagMax-lagGrid, 1+lagMax, 1+lagMax+lagGrid, 1+2*lagMax];
xcmapXticklabel = round( (xcmapXtick-1-lagMax)*MDtimeInterval_, 2);

fchanZ = cell(1, maxLayer);

for indL = 1:maxLayer

    mat2 = cell2mat(ch1Zsamples{indL});
    inputmap = mat2;
    %alphadata0 = ones(size(inputmap));
    %alphadata0(~Mask) = 0.2;
    %alphadata0(isnan(inputmap)) = 0;

    % nan case
    if all(isnan(inputmap))
        smmap = nan(size(inputmap));
    else
        smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
    end

    fchanZ{indL} = figure('Visible', p.figFlag);
    figtmp = imagesc(smmap);
    title([chanTitle, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)

    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;

    set(gca, 'XTick', xcmapXtick)
    set(gca, 'XTickLabel', {xcmapXticklabel})
    set(gca, 'XGrid', 'on')

    ax.FontSize = 15;

end

%%
for indL = 1:maxLayer
    saveas(fchanZ{indL}, fullfile(figuresDir, [chanName, '_minVel_ch1Zsamples', num2str(indL), 'L.png']), 'png')
    saveas(fchanZ{indL}, fullfile(figuresDir, [chanName, '_minVel_ch1Zsamples', num2str(indL), 'L.fig']), 'fig')
end

%% close figures
pause(0.5)
for indL = 1:maxLayer
    close(fchanZ{indL})
end
 
%%
disp('==== phaseDescriptives_MaxMinVel_OneChan is Done! ====')

   
end
