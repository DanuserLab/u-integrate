function [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, varargin)
% phaseMasking COMPUTE protrusion/retraction masks for an edge velocity map
% through over-smoothing and simple thresholding.
%
% Usage:
% [protMask1, retMask1] = phaseMasking(MD, smParamTh, figuresDir, varargin)
%
% Jungsik Noh, 2016/11/01
% Update by Qiongjing (Jenny) Zou, Aug 2018 
%
% Copyright (C) 2024, Danuser Lab - UTSouthwestern 
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
ip.addParameter('figFlag', 'off');
ip.addParameter('minimumRunLength', 5);    % Alternatively, you can remove 
                                           % index s.t. RunLength < 10
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);             
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 


parse(ip, varargin{:});
p = ip.Results;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end


%%  getting Maps from channels

chan0Title = 'Velocity (nm/sec)';
chan0Name = 'Vel';
disp(chan0Name)
disp(chan0Title)

iChan = 0;
maxLayer = 1;


[~, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
        'omittedWindows', p.omittedWindows, 'Folding', p.Folding, ...
        'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing, ...
        'EWMA', p.EWMA); 
    
disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])
% ..st layer

velmap = rawActmap{1}; size(velmap)
velmap_outl = actmap_outl{1}; size(velmap_outl)
imvelocitymap = imActmap{1}(:, 2:tmax);     %  Imputation (used as an input of computations)
                                                % Note 2:tmax


%%  smoothActivityMap prot/act maps

% smParamTh = 0.01;

inputmap = [nan(wmax, 1), imvelocitymap];

filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParamTh, 'UpSample', 1);
fsmParamTh = figure('Visible', p.figFlag);
figtmp = imagesc(filteredmap);
title([chan0Title, ' smParam=', num2str(smParamTh)])
colorbar;colormap(jet)

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%%
cmin = min(filteredmap(:)); disp(cmin)
cmax = max(filteredmap(:)); disp(cmax)

%% signF

signF = sign(filteredmap);

%%

fsignF = figure('Visible', p.figFlag);
 
inputmap = signF;
figtmp = imagesc(inputmap);
title('Sign of smoothedMap')
colorbar; 
cmap0 = jet;
colormap([cmap0(1,:); cmap0(end,:)])

figtmp.AlphaData = 1-isnan(inputmap);
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%%

protMask0 = (signF == 1);
retMask0 = (signF == -1);

%%
%protOnly = [nan(wmax, 1), imvelocitymap];
%protOnly(retMask0) = nan;


%%  small runs (<p.minimumRunLength) deletion
% minimumRunLength = 1;

wmax = size(protMask0, 1);
protMask1 = zeros(size(protMask0));

for w = 1:wmax
    xx = protMask0(w, :);
    xxrle = rle(xx);
    val = xxrle(1:2:end);
    vallength = xxrle(2:2:end);
    % threshold <- less than 10
    val(vallength < p.minimumRunLength) = 0;   

    xxrle2 = nan(size(xxrle));
    xxrle2(1:2:end) = val;
    xxrle2(2:2:end) = vallength;

    xxDel = irle(xxrle2);
    %tmp = [xx; xxDel];
    protMask1(w, :) = xxDel;
end


wmax = size(retMask0, 1);
retMask1 = zeros(size(retMask0));

for w = 1:wmax
    xx = retMask0(w, :);
    xxrle = rle(xx);
    val = xxrle(1:2:end);
    vallength = xxrle(2:2:end);
    % threshold <- less than 10
    val(vallength < p.minimumRunLength) = 0;   

    xxrle2 = nan(size(xxrle));
    xxrle2(1:2:end) = val;
    xxrle2(2:2:end) = vallength;

    xxDel = irle(xxrle2);
    %tmp = [xx; xxDel];
    retMask1(w, :) = xxDel;
end




%%  smoothActivityMap prot/act maps

smParam = smParamTh;

inputmap = [nan(wmax, 1), imvelocitymap];
alphadata0 = ones(size(inputmap));
alphadata0(~protMask1) = 0.2;

smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);

fprotSm = figure('Visible', p.figFlag);
figtmp = imagesc(smmap, [cmin, cmax]);
title('Protrution phase')
colorbar;colormap(jet)

%figtmp.AlphaData = 1-isnan(inputmap);
figtmp.AlphaData = alphadata0;
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%
smParam = 0.9;

inputmap = [nan(wmax, 1), imvelocitymap];
alphadata0 = ones(size(inputmap));
alphadata0(~protMask1) = 0.2;

smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
fprot = figure('Visible', p.figFlag);
figtmp = imagesc(smmap);
title('Protrution phase')
colorbar;colormap(jet)

%figtmp.AlphaData = 1-isnan(inputmap);
figtmp.AlphaData = alphadata0;
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%%  smoothActivityMap prot/act maps

smParam = smParamTh;

inputmap = [nan(wmax, 1), imvelocitymap];
alphadata0 = ones(size(inputmap));
alphadata0(~retMask1) = 0.2;

smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
fretSm = figure('Visible', p.figFlag);
figtmp = imagesc(smmap, [cmin, cmax]);
title('Retraction phase')
colorbar;colormap(jet)

%figtmp.AlphaData = 1-isnan(inputmap);
figtmp.AlphaData = alphadata0;
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;


%
smParam = 0.9;

inputmap = [nan(wmax, 1), imvelocitymap];
alphadata0 = ones(size(inputmap));
alphadata0(~retMask1) = 0.2;

smmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
fret = figure('Visible', p.figFlag);
figtmp = imagesc(smmap);
title('Retraction phase')
colorbar;colormap(jet)

%figtmp.AlphaData = 1-isnan(inputmap);
figtmp.AlphaData = alphadata0;
axis xy;xlabel('Time (s)');ylabel('Window')
ax = gca;
curTick = ax.XTick;
ax.XTickMode = 'manual';
ax.XTick = curTick+1;
ax.XTickLabel = (curTick)*MDtimeInterval_;

%%
saveas(fsmParamTh, fullfile(figuresDir, ['/fsmParamTh_', num2str(smParamTh), '.png']), 'png')
saveas(fsignF, fullfile(figuresDir, ['/fsignF_', num2str(smParamTh), '.png']), 'png')
saveas(fprotSm, fullfile(figuresDir, ['/fprotSm_', num2str(smParamTh), '.png']), 'png')
saveas(fprot, fullfile(figuresDir, ['/fprot_', num2str(smParamTh), '.png']), 'png')
saveas(fretSm, fullfile(figuresDir, ['/fretSm_', num2str(smParamTh), '.png']), 'png')
saveas(fret, fullfile(figuresDir, ['/fret_', num2str(smParamTh), '.png']), 'png')

%%
saveas(fsmParamTh, fullfile(figuresDir, ['/fsmParamTh_', num2str(smParamTh), '.fig']), 'fig')
saveas(fsignF, fullfile(figuresDir, ['/fsignF_', num2str(smParamTh), '.fig']), 'fig')
saveas(fprotSm, fullfile(figuresDir, ['/fprotSm_', num2str(smParamTh), '.fig']), 'fig')
saveas(fprot, fullfile(figuresDir, ['/fprot_', num2str(smParamTh), '.fig']), 'fig')
saveas(fretSm, fullfile(figuresDir, ['/fretSm_', num2str(smParamTh), '.fig']), 'fig')
saveas(fret, fullfile(figuresDir, ['/fret_', num2str(smParamTh), '.fig']), 'fig')


%%
disp('==== phaseMasking is Done! ====')

end
