function [protMask1, retMask1, filteredmap] = phaseMaskingInternal(MD, smParamTh, varargin)
% phaseMaskingInternal
% [protMask1, retMask1] = phaseMaskingInternal(MD, smParamTh, varargin)
%
% Updates:
% J Noh, 2019/07/05. Add 'EWMA' option. 
% J Noh, 2019/04/18. Add 'movingAvgSmoothing' option here too  to take 
% smoothed activity map input. 
% Jungsik Noh, 2017/03/30
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
ip.addParameter('minimumRunLength', 1);    % Alternatively, you can remove 
                                           % index s.t. RunLength < 10

ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);      
ip.addParameter('movingAvgSmoothing', false);   % 2019/04/18. For smoothed input
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing. 

parse(ip, varargin{:});
p = ip.Results;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
%if ~isdir(figuresDir); mkdir(figuresDir); end


%%  getting Maps from channels

chan0Title = 'Velocity (nm/sec)';
chan0Name = 'Vel';
disp(chan0Name)
disp(chan0Title)

iChan = 0;
maxLayer = 1;


[~, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
        'omittedWindows', p.omittedWindows, 'Folding', p.Folding, 'subFrames', p.subFrames, ...
        'movingAvgSmoothing', p.movingAvgSmoothing, 'EWMA', p.EWMA);
    
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

%filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParamTh, 'UpSample', 1);
% nan case
if all(isnan(inputmap))
    filteredmap = nan(size(inputmap));
else
    filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParamTh, 'UpSample', 1);
end


%%
%cmin = min(filteredmap(:)); disp(cmin)
%cmax = max(filteredmap(:)); disp(cmax)

%% signF

signF = sign(filteredmap);

%%

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


%%
disp('==== phaseMaskingInternal is Done! ====')

end
