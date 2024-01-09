function MD_iGC_SPAR2ch(MD, iChan1, iChan2,  chan1Name, chan2Name, layerMax, ...
    figuresDir, twlagMax, twlagMaxReg, varargin)
% MD_iGC_SPAR2ch() COMPUTE P-values for instantaneous Granger-Causality (iGC) 
% from Chan1 to Chan2 in each window. 
%
% Regression Model:
%   iChan2_t ~ lagged iChan2_t + lagged iChan1_t (AR+Reg)
%
% Usage example:
%   MD_iGC_SPAR2ch(MD, iChan1, iChan2, chan1Name, chan2Name, layerMax, ...
%    figuresDir, twlagMax, twlagMaxReg, 'WithN', 0, 'omittedWindows', omitWin, ...
%    'subFrames', subFr, 'parpoolNum', 10, ...
%    'movMedFrameSize', GCparam.movMedFrameSize, 'movingAvgSmoothing', GCparam.movingAvgSmoothing, ...
%    'CommonFactorNormAddChVec', p.CommonFactorNormAddChVec, 'factoranMethod', GCparam.factoranMethod, ...
%    'baseOfRatioVec', p.baseOfRatioVec, 'EWMA', GCparam.EWMA, ...
%    'infoCriterion', GCparam.infoCriterion)
%
% Updates:
% J Noh, 2021/07/21. Due to Java memory error in GUI, add figure close
% commands.
% J Noh, 2021/02/11. Rename and add comments. Output filename for
% '*_winFits_Layers.mat' is adjusted.
% J Noh,2020/06/16. Add .fig output for '_PvalPairedMaps_Vslztn' plot.
% J Noh, 2019/12/26. Collect partial-R-square output. 
% J Noh, 2019/10/01. Add smoothed figPairedAct2 for visualization purpose.
% J Noh, 2019/05/04. 'EWMA' specifies lambda value for exponentially
% weighted moving average smoothing. Default is 1 (raw data).
% J Noh, 2019/04/30. Add 'ratio' chan (4xx) and a 'baseOfRatio' option.
% J Noh, 2019/04/12. Can select BIC option via 'infoCriterion'.
% J Noh, 2019/04/02. 4-way: Left/Right + Up+Down propagating AR.
% J Noh, 2019/02/28. Add an option, 'factoranMethod'.
% J Noh, 2018/10/30. Remove the option 'CFnormalize' in calling
% mapOutlierImputation().
% J Noh, 2018/05/30. Add Common Factor (CF) analysis-based normalization.
% modified from GC_2chAdjAR.m 2018/05/30
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


tlagMax = twlagMax(1);
wlagMax = twlagMax(2);
tlagMaxReg = twlagMaxReg(1);
wlagMaxReg = twlagMaxReg(2);        % obsolete

% wlagMaxReg should be <= wlagMax
if (wlagMaxReg > wlagMax); wlagMaxReg = wlagMax; end

pcorrLagMax0 = floor(MD.nFrames_/10);

ip = inputParser;

ip.addParameter('figFlag', 'off');
ip.addParameter('impute', true);
ip.addParameter('WithN', false);
ip.addParameter('parpoolNum', 4);

ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('movingAvgSmoothing', false);
ip.addParameter('partialCorrLagMax', pcorrLagMax0)
ip.addParameter('movMedFrameSize', nan);   % movmedian frame size for normalization.
%ip.addParameter('CFnormalizeVec', [false, false]);
ip.addParameter('CommonFactorNormAddChVec', {NaN, NaN, NaN});
ip.addParameter('factoranMethod', 33);
ip.addParameter('infoCriterion', 'AIC');
ip.addParameter('baseOfRatioVec', [NaN, NaN, NaN]);  % base chan for ratio comp. default is ch1.
ip.addParameter('EWMA', 1);         % lambda=1 means no smoothing.


ip.parse(varargin{:});
p = ip.Results;

%p.figFlag = p.figFlag;
set(groot,'defaultLegendAutoUpdate','off')

%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isfolder(figuresDir); mkdir(figuresDir); end

tmptext = ['MD_GC_4wSPAR3ch_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')


%%  getting Maps from channels 1, 2

disp(chan1Name)


[~, ~,MDtimeInterval_, wmax, tmax, ~, ~, imActmap1] ...
    = mapOutlierImputation(MD, iChan1, layerMax+1, 'impute', p.impute, ...
    'omittedWindows', p.omittedWindows, 'WithN', p.WithN, 'subFrames', p.subFrames, ...
    'movingAvgSmoothing', p.movingAvgSmoothing, 'movMedFrameSize', p.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddChVec{1}, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatioVec(1), 'EWMA', p.EWMA, ...
    'figuresDir', figuresDir, 'chanName', chan1Name);

disp(chan2Name)

[~, ~, ~, ~, ~, ~, ~, imActmap2] ...
    = mapOutlierImputation(MD, iChan2, layerMax+1, 'impute', p.impute, ...
    'omittedWindows', p.omittedWindows, 'WithN', p.WithN, 'subFrames', p.subFrames, ...
    'movingAvgSmoothing', p.movingAvgSmoothing, 'movMedFrameSize', p.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddChVec{2}, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatioVec(2), 'EWMA', p.EWMA, ...
    'figuresDir', figuresDir, 'chanName', chan2Name);

fsaveName0 = ['frCh', num2str(iChan1), 'toCh', num2str(iChan2)];
disp(fsaveName0)
% 2021/02/11
fsaveName1 = ['fr_', chan1Name, '_to_', chan2Name];
disp(fsaveName1)

%%  to handle vel channel reads

if layerMax > 1
    if numel(imActmap1) == 1
        for indL = 2:layerMax; imActmap1{indL} = imActmap1{1}; end
    end
    if numel(imActmap2) == 1
        for indL = 2:layerMax; imActmap2{indL} = imActmap2{1}; end
    end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% due to parfor
%if isempty(gcp('nocreate')); parpool('local', p.parpoolNum); end


%%  global nanK, nanK2

nanK = 0;

for indL = 1:layerMax
    map = imActmap1{indL};
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK = max(nanK, map_ik);
end

for indL = 1:layerMax
    map = imActmap2{indL};
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK = max(nanK, map_ik);
end


disp(['nanK: ', num2str(nanK)])

%

nanK2 = 0;

for indL = 1:layerMax
    map = flip(imActmap1{indL}, 2);
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK2 = max(nanK2, map_ik);
end

for indL = 1:layerMax
    map = flip(imActmap2{indL}, 2);
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK2 = max(nanK2, map_ik);
end


disp(['nanK2: ', num2str(nanK2)])


%%  AR-Reg model selection for chan2 ~ lagged(chan2) + lagged(chan1) + given other chan(s)
% SPAR_4w3ch_cvLsoAIC, left/right/up/down, vel chan needs special handling.
% 201904, Parse MD.maps into SPAR-regression inputs
% 201904, Now drop out w1, wEnd connection for round cell

%tic
winch2Pvec = cell(layerMax, 1);    % AR orders
winch1Pvec = cell(layerMax, 1);    % orders for Regression terms

modelIndVarsMat = cell(layerMax, 1);
modelCoefMat = cell(layerMax, 1);

fmdlSlctn = cell(layerMax, 1);
Jw0s = cell(layerMax, 1); Jw1s = cell(layerMax, 1);
Pw0s = cell(layerMax, 1); Pw1s = cell(layerMax, 1);
Qw0s = cell(layerMax, 1);
resiMap = cell(layerMax, 1);
yhatFMap = cell(layerMax, 1);
yhatRMap = cell(layerMax, 1);
winFits = cell(layerMax, 1);

% AutoRegressive models with spatial propagation
% combination with repetition
%wlagMax_tmp = 1;

% for reg_coef_structure
%tmp = nmultichoosek(0:tlagMaxReg, wlagMaxReg+1);
% wlagMax===1
% UDLRAllLag01ctrledGC
tmp = nmultichoosek(1:tlagMax, 1+1);
ch2PvecMat = flip(tmp, 2);
ch1PvecMat = ch2PvecMat(1:end, :);      % except [0,0,0] = [Pwin0, Pwin1, Pwin2] to always include X_t.
% except [0,0,0] = [Pwin0, Pwin1, Pwin2] to always include X_t.

for indL = 1:layerMax
    
    % prepare maps, pre-Parsing
    % only adjacent layers (up, down)
    mapy = cell(3,1); mapx = cell(3,1); %mapz = cell(3,1);
    
    mapy{2} = imActmap2{indL};
    mapy{1} = nan(size(mapy{2})); mapy{3} = nan(size(mapy{2}));
    if (iChan2 ~= 0) && (indL > 1); mapy{1} = imActmap2{indL-1}; end    % up
    if (iChan2 ~= 0); mapy{3} = imActmap2{indL+1}; end                  % down
    
    mapx{2} = imActmap1{indL};
    mapx{1} = nan(size(mapx{2})); mapx{3} = nan(size(mapx{2}));
    if (iChan1 ~= 0) && (indL > 1); mapx{1} = imActmap1{indL-1}; end    % up
    if (iChan1 ~= 0); mapx{3} = imActmap1{indL+1}; end                  % down
    
%    mapz{2} = imActmap3{indL};
%    mapz{1} = nan(size(mapz{2})); mapz{3} = nan(size(mapz{2}));
%    if (iChan3 ~= 0) && (indL > 1); mapz{1} = imActmap3{indL-1}; end    % up
%    if (iChan3 ~= 0); mapz{3} = imActmap3{indL+1}; end                  % down
    
    
    %mapy_ik = find(~all(isnan(mapy)), 1) - 1;
    %if isempty(mapy_ik); mapy_ik = 0; end   % for the case of map is all nans.
    %mapx_ik = find(~all(isnan(mapx)), 1) - 1;
    %if isempty(mapx_ik); mapx_ik = 0; end   % for the case of map is all nans.
    %nanK = max(mapx_ik, mapy_ik);
    
    if  (nanK > 0) || (nanK2 > 0)
        for k = 1:3
            mapx{k} = mapx{k}(:, nanK+1:tmax-nanK2);
            mapy{k} = mapy{k}(:, nanK+1:tmax-nanK2);
             
        end
    end
    
    % standardization for lasso
    zmapy = cell(3,1); zmapx = cell(3,1); %zmapz = cell(3,1);
    for k = 1:3
        zmapx{k} = zscore(mapx{k}')';
        zmapy{k} = zscore(mapy{k}')'; 
    end
    
    % adj AR-Reg for all models
    % it does Parsing + LsoReg + AICselection, inputmap has 3 layers
    [arMdlAvBICmap, arMdlAvBICcurve, reducedMdlAvBICmap, reducedMdlAvBICcurve, ...
        fullMdlAvBICmap, fullMdlAvBICcurve, Jw0, Jw1, ...
        Pw1, Qw0, modelIndVars_mat, modelCoef_mat, winFits_mat, resi_mat, yhatF_mat, yhatR_mat] = ...
    SPAR_2ch_mdlSlctn_iGC_3steps(zmapx, zmapy, tlagMax, wlagMax, ...
        tlagMaxReg, wlagMaxReg, ch1PvecMat, p.infoCriterion);
    
    fmdlSlctn{indL} = figure('Visible', p.figFlag);
    subplot(3,2,1)
    plot(arMdlAvBICmap')
    title([num2str(indL), 'L ', 'Per-window ICs'])
    xlabel('Model orders')
    ylabel('Avg IC')
    
    subplot(3,2,2)
    plot(arMdlAvBICcurve)
    title([num2str(indL), 'L ', 'AR Mdl: Jw0 Jw1= ', num2str(Jw0), ', ', num2str(Jw1)])
    xlabel('Model orders')
    ylabel('Avg IC')  
    
    subplot(3,2,3)
    plot(reducedMdlAvBICmap')
    title([num2str(indL), 'L ', 'Per-window ICs'])
    xlabel('Model orders')
    ylabel('Avg IC')
    
    subplot(3,2,4)
    plot(reducedMdlAvBICcurve)
    title([num2str(indL), 'L ', 'Ctrl Mdl: Pw1= ', num2str(Pw1)])
    xlabel('Model orders')
    ylabel('Avg IC')    
    
    subplot(3,2,5)
    plot(fullMdlAvBICmap')
    title([num2str(indL), 'L ', 'Per-window ICs'])
    xlabel('Model orders')
    ylabel('Avg IC')
    
    subplot(3,2,6)
    plot(fullMdlAvBICcurve)
    title([num2str(indL), 'L ', 'Full Mdl: Qw0= ', num2str(Qw0)])
    xlabel('Model orders')
    ylabel('Avg IC')       
    
    saveas(fmdlSlctn{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_fmdlSlctn.png']), 'png')
    saveas(fmdlSlctn{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_fmdlSlctn.fig']), 'fig')
    %end
    
    modelIndVarsMat{indL} = modelIndVars_mat;
    modelCoefMat{indL} = modelCoef_mat;
    Jw0s{indL} = Jw0;
    Jw1s{indL} = Jw1;
    %Pw0s{indL} = Pw0;
    Pw1s{indL} = Pw1;
    Qw0s{indL} = Qw0;
    resiMap{indL} = resi_mat;
    yhatFMap{indL} = yhatF_mat;
    yhatRMap{indL} = yhatR_mat;
    winFits{indL} = winFits_mat;
    %armaxLagMat{indL} = tmpmaxLagsMat;
    %winch2Pvec{indL} = ch2Pvectmp;
    %winch1Pvec{indL} = ch1Pvectmp;
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(fmdlSlctn{indL})
end

%%  plot model selection by Lasso + OLS AIC

figwinIndMat = cell(layerMax, 1);
%xtickvec = [(1:2*(1+2*wlagMax+2)) * tlagMax, 2*(1+2*wlagMax+2)*tlagMax + (1:(1+2*wlagMaxReg+2)) * tlagMaxReg];
%xtickvec = [(1:1*(1+2*wlagMax+2)) * tlagMax, (6:2*(1+2*wlagMax+2))*tlagMax + 1, ...
%    (1+2*(1+2*wlagMax+2))*tlagMax+1:tlagMaxReg:2*(1+2*wlagMax+2)*tlagMax+1+tlagMaxReg*(2+wlagMaxReg*2), ...
%    2*(1+2*wlagMax+2)*tlagMax+1+tlagMaxReg*(1+2+wlagMaxReg*2)+1];
%xtickvec = [(1:5)*tlagMax, 5*tlagMax + (1:10)*(tlagMax+1)];
xtickvec = [(1:5)*tlagMax, 5*tlagMax + (1:5)*(tlagMax+1)];


for indL = 1:layerMax
    
    figwinIndMat{indL} = figure('Visible', p.figFlag);
    figtmp = imagesc(modelIndVarsMat{indL});
    figtmp.AlphaData = 1-isnan(modelIndVarsMat{indL});
    %
    axis xy
    ax = gca;
    ax.XTick = xtickvec;
    avgArmaxLag = round(mean(winch2Pvec{indL}, 'omitnan'), 0);
    avgRegLag = round(mean(winch1Pvec{indL}, 'omitnan'), 0);
    
    avgMaxOrder = [Jw0s{indL}, Jw1s{indL}, Pw0s{indL}, Pw1s{indL}, Qw0s{indL}];
    
    hold on
    for i=xtickvec
        line([i+0.5,i+0.5], [1,wmax])
    end
    
    title({['AR/Ctrl/Reg features chosen by IC - '];[ fsaveName0, '-', num2str(indL), 'L']; ...
        ['Orders : ', sprintf('%.0f ',avgMaxOrder)]})
    ylabel('Windows')
    
    %legend('Pwin0', 'Pwin+1', 'Pwin+2', 'Pwin+3', 'Location','northoutside','Orientation','horizontal');
    
    ax = gca; ax.FontSize = 10;
    
    %%
    saveas(figwinIndMat{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_modelIndVarsMat.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figwinIndMat{indL})
end

%%   modelCoefMat

figwinCoefMat = cell(layerMax, 1);
%xtickvec = [(1:2*(1+2*wlagMax+2)) * tlagMax, 2*(1+2*wlagMax+2)*tlagMax + (1:(1+2*wlagMaxReg+2)) * tlagMaxReg];

for indL = 1:layerMax
    
    figwinCoefMat{indL} = figure('Visible', p.figFlag);
    
    inputmap = atan(modelCoefMat{indL});
    %inputmap(inputmap == 0) = NaN;
    figtmp = imagesc(inputmap);
    figtmp.AlphaData = 1-isnan(inputmap);
    colormap(jet); colorbar
    maxval = max(abs(inputmap(:)));
    if ~isnan(maxval); caxis([-maxval, maxval]); end
    %
    axis xy
    ax = gca;
    ax.XTick = xtickvec;
    avgArmaxLag = round(mean(winch2Pvec{indL}, 'omitnan'), 0);
    avgRegLag = round(mean(winch1Pvec{indL}, 'omitnan'), 0);
    
    avgMaxOrder = [Jw0s{indL}, Jw1s{indL}, Pw0s{indL}, Pw1s{indL}, Qw0s{indL}];
    
    hold on
    for i=xtickvec
        line([i+0.5,i+0.5], [1,wmax])
    end
    
    title({['AR/Ctrl/Reg features chosen by IC - '];[ fsaveName0, '-', num2str(indL), 'L']; ...
        ['Order: ', sprintf('%.0f ',avgMaxOrder)]})
    ylabel('Windows')
    
    %legend('Pwin0', 'Pwin+1', 'Pwin+2', 'Pwin+3', 'Location','northoutside','Orientation','horizontal');
    
    ax = gca; ax.FontSize = 10;
    
    %%
    saveas(figwinCoefMat{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_figwinCoefMat.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figwinCoefMat{indL})
end

%%
save(fullfile(figuresDir, [fsaveName0, '_modelIndVarsMat_Layers.mat']), 'modelIndVarsMat')
save(fullfile(figuresDir, [fsaveName0, '_modelCoefMat_Layers.mat']), 'modelCoefMat')
save(fullfile(figuresDir, [fsaveName0, '_chosenOrder_Layers.mat']), 'Pw0s', 'Pw1s', 'Qw0s')


%%  plot mean curve of modelCoef

figCoefcurve = cell(layerMax, 1);
%xtickvec = [(1:2*(1+2*wlagMax+2)) * tlagMax, ...
%    2*(1+2*wlagMax+2)*tlagMax + (1:(1+2*wlagMaxReg+2)) * tlagMaxReg];
tickMax = xtickvec(end);

for indL = 1:layerMax
    figCoefcurve{indL} = figure('Visible', p.figFlag);
    
    inputmap = modelCoefMat{indL};
    coefCurve = mean(inputmap, 1, 'omitnan');
    coefAR = coefCurve; %coefAR((1+2*wlagMax+2)*tlagMax + 1:tickMax) = NaN;
    coefAR(xtickvec(5)+1:xtickvec(end)) = NaN;
    bar(coefAR, 'k');
    hold on
    coefCtrl = coefCurve; %coefCtrl(1:(1+2*wlagMax+2)*tlagMax) = NaN;
    %coefCtrl(2*(1+2*wlagMax+2)*tlagMax + (2*wlagMaxReg+2)*tlagMaxReg + 1:tickMax) = NaN;
    coefCtrl(1:xtickvec(1+2*wlagMax+2)) = NaN;
    coefCtrl(xtickvec(end-1)+1:tickMax) = NaN;
    bar(coefCtrl, 'b');
    coefReg = coefCurve; %coefReg(1:2*(1+2*wlagMax+2)*tlagMax + (2*wlagMaxReg+2)*tlagMaxReg) = NaN;
    coefReg(1:xtickvec(end-1)) = NaN;
    bar(coefReg, 'r');
    
    %    bar(coefCurve)
    h = refline([0,0]); h.Color = [.3 .3 .3];
    set(gca, 'XTick', xtickvec);
    for i = xtickvec
        l = line([i+0.5,i+0.5], ylim); l.Color = [.3 .3 .3];
    end
    
    avgArmaxLag = round(mean(winch2Pvec{indL}, 'omitnan'), 0);
    avgRegLag = round(mean(winch1Pvec{indL}, 'omitnan'), 0);
    avgMaxOrder = [Jw0s{indL}, Jw1s{indL}, Pw0s{indL}, Pw1s{indL}, Qw0s{indL}];
    
    title({['AR/Ctrl/Reg features chosen by IC - '];[ fsaveName0, '-', num2str(indL), 'L']; ...
        ['Order: ', sprintf('%.0f ',avgMaxOrder)]})
    ylabel('Mean coef')
    %legend('Pwin0', 'Pwin+1', 'Pwin+2', 'Pwin+3', 'Location','northoutside','Orientation','horizontal');
    ax = gca; ax.FontSize = 10;
    
    %
    saveas(figCoefcurve{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_figCoefcurve.png']), 'png')
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figCoefcurve{indL})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Partial corr map/curve in AIC-chosen adjAR-Reg models  2019/04 (deleted)


%%  plot orders by AIC: AR, chan2
%%  plot orders by AIC: Reg, chan1

%% Chosen fit, Ftest, FPval, Rsq, cnum, sigma, resiMap
% pval, partialCorr, resiMap, Rsq, cnum, normality, ~, ~
%% chosen fit
%% saveas winFits (2021/02/11 output file name changed)

save(fullfile(figuresDir, ['GC_', fsaveName1, '_winFits_Layers.mat']), 'winFits', ...
    'resiMap','yhatFMap','yhatRMap')


%% draw outcome
%% draw ECDF of Fpval, LRTpval, normalityTest


figecdf = cell(layerMax, 1);

for indL = 1:layerMax
    figecdf{indL} = figure('Visible', p.figFlag);
    
    Fpvec = winFits{indL}(:, 8);
    LRTpvec = winFits{indL}(:, 10);
    KSpvec = winFits{indL}(:, 5);
    LBpvec = winFits{indL}(:, 12);
    
    if ~all(isnan(Fpvec))
        [f,x] = ecdf(Fpvec);
        stairs(x, f)%, 'LineWidth', 2)
        
        
        hold on
        
        [f,x] = ecdf(LRTpvec);
        stairs(x, f, ':')%, 'LineWidth', 2)
        [f,x] = ecdf(KSpvec);
        stairs(x, f)
        [f,x] = ecdf(LBpvec);
        stairs(x, f)
    else
        x = nan; f=nan; stairs(x,f)
    end
    
    legend('F-test', 'LRT', 'resiNormality', 'LB-test', 'Location','southeast','Orientation','vertical');
    
    h=refline([1,0]); %h.Color='k';
    h=line([0.05, 0.05], [0,1]); h.Color='red';
    %h=refline([0, 0.5]);h.Color=[0.5 0.5 0.5];
    
    winProp = sum(Fpvec < 0.05)/sum(~isnan(Fpvec));
    medFPval = median(Fpvec, 'omitnan');
    
    subtitle0 = ['medianPval: ', num2str(round(medFPval, 3)), ', prop of <5%: ', num2str(round(winProp,2))];
    title0 = ['GC P-values from ', chan1Name, ' to ', chan2Name, ' -', num2str(indL), 'L'];
    title({title0; subtitle0})
    
    xlabel('P-value')
    ylabel('Cumulative Relative Frequency')
    
    ax = gca; ax.FontSize = 13;
    
    % save
    saveas(figecdf{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_GCecdf.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figecdf{indL})
end

%% paired Activity maps

%smParam = 0.8;      % 0.9 or 0.8
smParam = 1;        % since 2019/05/04. After SNR

figPairedAct = cell(layerMax, 1);

for indL = 1:layerMax
    
    Xmap = imActmap1{indL};
    Ymap = imActmap2{indL};
    
    if ~all(isnan(Xmap(:))) && ~all(isnan(Ymap(:)))
        
        smXmap = smoothActivityMap(Xmap, 'SmoothParam', smParam, 'UpSample', 1);
        smYmap = smoothActivityMap(Ymap, 'SmoothParam', smParam, 'UpSample', 1);
        
        figPairedAct{indL} = figure('Visible', p.figFlag, 'Position', [317 330 7*80*2 420]);
        
        subplot(1,3,1);
        figtmp = imagesc(smXmap, quantile(smXmap(:), [0.001, 0.999]));
        title([chan1Name, '-', num2str(indL), 'L'])
        colormap(jet); colorbar
        figtmp.AlphaData = 1-isnan(Xmap);
        axis xy;xlabel('Time (s)');ylabel('Window')
        ax = gca;
        curTick = ax.XTick;
        ax.XTickMode = 'manual';
        ax.XTick = curTick+1;
        ax.XTickLabel = (curTick)*MDtimeInterval_;
        
        %
        subplot(1,3,2);
        figtmp = imagesc(smYmap, quantile(smYmap(:), [0.001, 0.999]));
        title([chan2Name, '-', num2str(indL), 'L'])
        colormap(jet); colorbar
        figtmp.AlphaData = 1-isnan(Ymap);
        axis xy;xlabel('Time (s)');%ylabel('Window')
        ax = gca;
        curTick = ax.XTick;
        ax.XTickMode = 'manual';
        ax.XTick = curTick+1;
        ax.XTickLabel = (curTick)*MDtimeInterval_;
        
        %
        Fpvec = winFits{indL}(:, 8);
        subplot(1,3,3);
        sigGC = -log10(Fpvec);
        ind1 = (sigGC > -log10(0.05));
        sigGC1 = sigGC;
        sigGC1(~ind1) = NaN;
        sigGC0 = sigGC;
        sigGC0(ind1) = NaN;
        
        b = barh(1:wmax, sigGC0, 'b');
        ax = gca; ax.YLim = [1, wmax];
        axis xy
        hold on
        b = barh(1:wmax, sigGC1, 'r');
        ax = gca; ax.YLim = [1, wmax];
        
        h = line([-log10(0.05), -log10(0.05)], [1, wmax]);
        ptick = [0.05, 0.01, 0.001, 0.0001];
        logptick = -log10(ptick);
        ax = gca;
        ax.XLim = [0, 5];
        ax.XTickMode = 'manual';
        ax.XTick = logptick;
        ax.XTickLabel = ptick;    %{'0.05', '0.01', '0.001', ' '};
        xlabel('P-value (-log10)')
        
        ax = gca; ax.FontSize = 13;
        
    else
        figPairedAct{indL} = figure('Visible', p.figFlag, 'Position', [317 330 7*80*2 420]);
    end
end

%% saveas
for indL = 1:layerMax
    saveas(figPairedAct{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_PvalPairedMaps.png']), 'png')
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figPairedAct{indL})
end

%% paired Activity maps

smParam = 0.8;      % 0.9 or 0.8
%smParam = 1;        % since 2019/05/04. After SNR

figPairedAct2 = cell(layerMax, 1);

for indL = 1:layerMax
    
    Xmap = imActmap1{indL};
    Ymap = imActmap2{indL};
    
    if ~all(isnan(Xmap(:))) && ~all(isnan(Ymap(:)))
        
        smXmap = smoothActivityMap(Xmap, 'SmoothParam', smParam, 'UpSample', 1);
        smYmap = smoothActivityMap(Ymap, 'SmoothParam', smParam, 'UpSample', 1);
        
        figPairedAct2{indL} = figure('Visible', p.figFlag, 'Position', [317 330 7*80*2 420]);
        
        subplot(1,3,1);
        figtmp = imagesc(smXmap, quantile(smXmap(:), [0.001, 0.999]));
        title([chan1Name, '-', num2str(indL), 'L'])
        colormap(jet); colorbar
        figtmp.AlphaData = 1-isnan(Xmap);
        axis xy;xlabel('Time (s)');ylabel('Window')
        ax = gca;
        curTick = ax.XTick;
        ax.XTickMode = 'manual';
        ax.XTick = curTick+1;
        ax.XTickLabel = (curTick)*MDtimeInterval_;
        
        %
        subplot(1,3,2);
        figtmp = imagesc(smYmap, quantile(smYmap(:), [0.001, 0.999]));
        title([chan2Name, '-', num2str(indL), 'L'])
        colormap(jet); colorbar
        figtmp.AlphaData = 1-isnan(Ymap);
        axis xy;xlabel('Time (s)');%ylabel('Window')
        ax = gca;
        curTick = ax.XTick;
        ax.XTickMode = 'manual';
        ax.XTick = curTick+1;
        ax.XTickLabel = (curTick)*MDtimeInterval_;
        
        %
        Fpvec = winFits{indL}(:, 8);
        subplot(1,3,3);
        sigGC = -log10(Fpvec);
        ind1 = (sigGC > -log10(0.05));
        sigGC1 = sigGC;
        sigGC1(~ind1) = NaN;
        sigGC0 = sigGC;
        sigGC0(ind1) = NaN;
        
        b = barh(1:wmax, sigGC0, 'b');
        ax = gca; ax.YLim = [1, wmax];
        axis xy
        hold on
        b = barh(1:wmax, sigGC1, 'r');
        ax = gca; ax.YLim = [1, wmax];
        % 2020/06/03
        %gridvals = -log10([1:-0.1:0.1, 0.1:-0.01:0.01, 0.01:-0.001:0.001, 0.001:-0.0001:0.0001, 0.0001:-0.00001:0.00001]);
        gridvals = -log10([1:-0.25:0.1, 0.1:-0.025:0.01, 0.01:-0.0025:0.001, 0.001:-0.00025:0.0001, 0.0001:-0.000025:0.00001]);
        for k = 1:numel(gridvals)
            h = line([gridvals(k), gridvals(k)], [1, wmax]); 
            h.Color = [.5 .5 .5, 0.4]; %h.LineStyle = '--';
        end        
                
        h = line([-log10(0.05), -log10(0.05)], [1, wmax]);
        ptick = [0.1, 0.05, 0.01, 0.001, 0.0001];
        logptick = -log10(ptick);
        ax = gca;
        ax.XLim = [0, 5];
        ax.XTickMode = 'manual';
        ax.XTick = logptick;
        ax.XTickLabel = ptick;    %{'0.05', '0.01', '0.001', ' '};
        xlabel('P-value')
        
        ax = gca; ax.FontSize = 13;
        % 2020/06/03
        ax.XTickLabelRotation = 45;
        
    else
        figPairedAct2{indL} = figure('Visible', p.figFlag, 'Position', [317 330 7*80*2 420]);
    end
end

%% saveas
for indL = 1:layerMax
    saveas(figPairedAct2{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_PvalPairedMaps_Vslztn.png']), 'png')
    saveas(figPairedAct2{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_PvalPairedMaps_Vslztn.fig']), 'fig')
    % some java error happened.
    %saveas(figPairedAct2{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_PvalPairedMaps_Vslztn.pdf']), 'pdf')
    print(figPairedAct2{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL), 'L_PvalPairedMaps_Vslztn.pdf']), ...
        '-dpdf', '-bestfit')
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(figPairedAct2{indL})
end

%% resiMap2 (adjusted before)
%%
save(fullfile(figuresDir, ['GC_', fsaveName0, '_resiMap_Layers.mat']), 'resiMap')
save(fullfile(figuresDir, ['GC_', fsaveName0, '_yhatFMap_Layers.mat']), 'yhatFMap')
save(fullfile(figuresDir, ['GC_', fsaveName0, '_yhatRMap_Layers.mat']), 'yhatRMap')


%% draw y yhat map

smParam = 1;

fyMap = cell(1, layerMax);
fyhatFMap = cell(1, layerMax);
fyhatRMap = cell(1, layerMax);

for indL = 1:layerMax
    
    tmp = imActmap2{indL}(:, nanK+1:tmax-nanK2);
    inputmap = [nan(wmax, nanK), zscore(tmp')', nan(wmax, nanK2)];
    
    switch all(isnan(inputmap(:)))
        case true
            filteredmap = nan(size(inputmap));
            fyMap{indL} = figure('Visible', p.figFlag);
            figtmp = imagesc(filteredmap);
        case false
            filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
            fyMap{indL} = figure('Visible', p.figFlag);
            %figtmp = imagesc(filteredmap);
            clim0 = quantile(abs(filteredmap(:)), 0.998);
            figtmp = imagesc(filteredmap, [-clim0, clim0]);
    end
    
    title(['Y of ', fsaveName0, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)
    
    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    ax = gca; ax.FontSize = 10;
    %
    inputmap = yhatFMap{indL};
    
    switch all(isnan(inputmap(:)))
        case true
            filteredmap = nan(size(inputmap));
            fyhatFMap{indL} = figure('Visible', p.figFlag);
            figtmp = imagesc(filteredmap);
        case false
            filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
            fyhatFMap{indL} = figure('Visible', p.figFlag);
            %figtmp = imagesc(filteredmap);
            figtmp = imagesc(filteredmap, [-clim0, clim0]);
    end
    
    title(['Yhat Fullmodel of ', fsaveName0, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)
    
    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    ax = gca; ax.FontSize = 10;
    %
    inputmap = yhatRMap{indL};
    
    switch all(isnan(inputmap(:)))
        case true
            filteredmap = nan(size(inputmap));
            fyhatRMap{indL} = figure('Visible', p.figFlag);
            figtmp = imagesc(filteredmap);
        case false
            filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
            fyhatRMap{indL} = figure('Visible', p.figFlag);
            %figtmp = imagesc(filteredmap);
            figtmp = imagesc(filteredmap, [-clim0, clim0]);
    end
    
    title(['Yhat Reducedmodel of ', fsaveName0, '-', num2str(indL), 'L'])
    colorbar;colormap(jet)
    
    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    ax = gca; ax.FontSize = 10;
    
    %
    saveas(fyMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yMap.png']), 'png')
    saveas(fyhatFMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yhatFMap.png']), 'png')
    saveas(fyhatRMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yhatRMap.png']), 'png')
    
    saveas(fyMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yMap.fig']), 'fig')
    saveas(fyhatFMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yhatFMap.fig']), 'fig')
    saveas(fyhatRMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_yhatRMap.fig']), 'fig')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(fyMap{indL})
    close(fyhatFMap{indL})
    close(fyhatRMap{indL})
end

%% draw resiMap, yhatF,R map

smParam = 1;

fresiMap = cell(1, layerMax);
for indL = 1:layerMax
    
    inputmap = resiMap{indL};
    
    switch all(isnan(inputmap(:)))
        case true
            filteredmap = nan(size(inputmap));
            fresiMap{indL} = figure('Visible', p.figFlag);
            figtmp = imagesc(filteredmap);
        case false
            filteredmap = smoothActivityMap(inputmap, 'SmoothParam', smParam, 'UpSample', 1);
            fresiMap{indL} = figure('Visible', p.figFlag);
            %figtmp = imagesc(filteredmap);
            figtmp = imagesc(filteredmap, quantile(filteredmap(:), [0.001, 0.999]));
    end
    
    title(['GC-Residuals of ', fsaveName0, '-', num2str(indL), 'L'])
    
    colorbar;colormap(jet)
    
    figtmp.AlphaData = 1-isnan(inputmap);
    axis xy;xlabel('Time (s)');ylabel('Window')
    ax = gca;
    curTick = ax.XTick;
    ax.XTickMode = 'manual';
    ax.XTick = curTick+1;
    ax.XTickLabel = (curTick)*MDtimeInterval_;
    
    ax = gca; ax.FontSize = 10;
    
    %%
    saveas(fresiMap{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_resiMap2.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(fresiMap{indL})
end

%%
%% draw autoCorr curves (means) acCurve = autoCorrCurvePermTest

numPerm = 200;
parpoolNum = p.parpoolNum;
rseed = 'shuffle';

acCurve = cell(1, layerMax);
for indL = 1:layerMax
    
    title0 = [fsaveName0, '-residuals-', num2str(indL), 'L'];
    [acCurve{indL}, ~] = autoCorrCurvePermTest_mean(resiMap{indL}, title0, MDtimeInterval_, ...
        numPerm, parpoolNum, rseed, 'figFlag', p.figFlag);
    
    saveas(acCurve{indL}, fullfile(figuresDir, ['acCurve_', title0, '.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(acCurve{indL})
end

%%  partial Rsquare

hist_rsq = cell(1, layerMax);

for indL = 1:layerMax
    
    Fpvec = winFits{indL}(:, 8);
    partR2 = winFits{indL}(:, 11);
    ind0 = (Fpvec >= 0.05);
    %
    partR2sig = partR2;
    partR2sig(ind0) = NaN;
    
    hist_rsq{indL} = figure('Visible', p.figFlag);
    %histogram(partR2, 'BinMethod', 'fd')
    plot(1:wmax, partR2)
    hold on
    scatter(1:wmax, partR2sig, [], 'r')
    refline([0, mean(partR2, 'omitnan')]);
    title(['partial-Rsquares: ', fsaveName0, '-', num2str(indL), 'L', ' mean: ', num2str(mean(partR2, 'omitnan'))])
    xlabel('R-squares')
    
    ax = gca; ax.FontSize = 10;
    
    %%
    saveas(hist_rsq{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_partialRsquares.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(hist_rsq{indL})
end

%%  condition number

hist_cnum = cell(1, layerMax);

for indL = 1:layerMax
    
    %Fpvec = winFits{indL}(:, 8);
    cnumVec = winFits{indL}(:, 2);
    %ind0 = (Fpvec >= 0.05);
    %partR2(ind0) = NaN;
    % imaginary number!?
    if ~isreal(cnumVec); cnumVec = nan(size(cnumVec)); end
    
    hist_cnum{indL} = figure('Visible', p.figFlag);
    %histogram(cnumVec, 'BinMethod', 'fd')
    plot(cnumVec)
    refline([0, mean(cnumVec, 'omitnan')]);
    h=refline([0, 30]); h.Color = 'r';
    title({['condNum fullmodAllwin: ', fsaveName0, '-', num2str(indL), 'L']; ...
        [' mean: ', num2str(mean(cnumVec, 'omitnan'))]})
    xlabel('condNum')
    
    ax = gca; ax.FontSize = 10;
    ax.YLim(1) = 0;
    
    %%
    saveas(hist_cnum{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_condNum.png']), 'png')
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(hist_cnum{indL})
end

%%
disp('====End of MD_iGC_SPAR2ch ====')

end

