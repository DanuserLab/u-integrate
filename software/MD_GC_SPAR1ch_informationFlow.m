function MD_GC_SPAR1ch_informationFlow(MD, iChan2, chan2Name, layerMax, ...
    figuresDir, twlagMax, varargin)
% MD_GC_SPAR1ch_informationFlow() COMPUTE P-values for intra-cellular information flows 
% based on Granger-Causality (GC) for one channel. 
%
% Updates:
% J Noj, 2021/07/21. Due to Java memory error in GUI, add figure close
% commands.
% J Noh, 2021/02/05. Rename and add comments. 
% Modified from MD_GC_SPAR2ch_UDLRAllLag0ctrledGC_allBIC3steps.m 
% J NOH. 2019/12/27.
%
% Copyright (C) 2021, Danuser Lab - UTSouthwestern 
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

disp(chan2Name)


[~, ~,MDtimeInterval_, wmax, tmax, ~, ~, imActmap2] ...
    = mapOutlierImputation(MD, iChan2, layerMax+1, 'impute', p.impute, ...
    'omittedWindows', p.omittedWindows, 'WithN', p.WithN, 'subFrames', p.subFrames, ...
    'movingAvgSmoothing', p.movingAvgSmoothing, 'movMedFrameSize', p.movMedFrameSize, ...
    'CommonFactorNormAddCh', p.CommonFactorNormAddChVec{1}, 'factoranMethod', p.factoranMethod, ...
    'baseOfRatio', p.baseOfRatioVec(1), 'EWMA', p.EWMA, ...
    'figuresDir', figuresDir, 'chanName', chan2Name);


fsaveName0 = ['intraCh', num2str(iChan2)];
disp(fsaveName0)

%%  to handle vel channel reads

if layerMax > 1
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
    map = imActmap2{indL};
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK = max(nanK, map_ik);
end

disp(['nanK: ', num2str(nanK)])

%

nanK2 = 0;

for indL = 1:layerMax
    map = flip(imActmap2{indL}, 2);
    map_ik = find(~all(isnan(map)), 1) - 1;
    if isempty(map_ik); map_ik = 0; end   % for the case of map is all nans.
    nanK2 = max(nanK2, map_ik);
end

disp(['nanK2: ', num2str(nanK2)])


%%  AR-Reg model selection for chan2 ~ lagged(chan2) (pure SPAR)
% SPAR_4w3ch_cvLsoAIC, left/right/up/down, vel chan needs special handling.
% 201904, Parse MD.maps into SPAR-regression inputs
% 201904, Now drop out w1, wEnd connection for round cell

%tic
winch2Pvec = cell(layerMax, 1);    % AR orders
%winch1Pvec = cell(layerMax, 1);    % orders for Regression terms

modelIndVarsMat = cell(layerMax, 1);
modelCoefMat = cell(layerMax, 1);

fmdlSlctn = cell(layerMax, 1);
Jw0s = cell(layerMax, 1); Jw1s = cell(layerMax, 1);
Pw0s = cell(layerMax, 1); Pw1s = cell(layerMax, 1);
Qw0s = cell(layerMax, 1);
resiMap = cell(layerMax, 1);
yhatFMap = cell(layerMax, 1);
yhatRMap = cell(layerMax, 1);
% 2019/12/30. LRUDgc
winFits = cell(1,4);
for direction = 1:4
    winFits{direction} = cell(layerMax, 1);
end

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
    
%    mapx{2} = imActmap1{indL};
%    mapx{1} = nan(size(mapx{2})); mapx{3} = nan(size(mapx{2}));
%    if (iChan1 ~= 0) && (indL > 1); mapx{1} = imActmap1{indL-1}; end    % up
%    if (iChan1 ~= 0); mapx{3} = imActmap1{indL+1}; end                  % down
    
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
            %mapx{k} = mapx{k}(:, nanK+1:tmax-nanK2);
            mapy{k} = mapy{k}(:, nanK+1:tmax-nanK2);
             
        end
    end
    
    % standardization for lasso
    zmapy = cell(3,1); zmapx = cell(3,1); %zmapz = cell(3,1);
    for k = 1:3
        %zmapx{k} = zscore(mapx{k}')';
        zmapy{k} = zscore(mapy{k}')'; 
    end
    
    % adj AR-Reg for all models
    % it does Parsing + LsoReg + AICselection, inputmap has 3 layers
    [arMdlAvBICmap, arMdlAvBICcurve,  ...
        Jw0, Jw1, ...
        modelIndVars_mat, modelCoef_mat, winFits_mat] = ...
    SPAR_1ch_mdlSlctn_UDLR_intraGC(zmapy, tlagMax, wlagMax, ...
         ch1PvecMat, p.infoCriterion);
    
    fmdlSlctn{indL} = figure('Visible', p.figFlag);
    subplot(1,2,1)
    plot(arMdlAvBICmap')
    title([num2str(indL), 'L ', 'Per-window ICs'])
    xlabel('Model orders')
    ylabel('Avg IC')
    
    subplot(1,2,2)
    plot(arMdlAvBICcurve)
    title([num2str(indL), 'L ', 'AR Mdl: Jw0 Jw1= ', num2str(Jw0), ', ', num2str(Jw1)])
    xlabel('Model orders')
    ylabel('Avg IC')  
    
   
    saveas(fmdlSlctn{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_fmdlSlctn.png']), 'png')
    saveas(fmdlSlctn{indL}, fullfile(figuresDir, [fsaveName0,'_',num2str(indL),'L_fmdlSlctn.fig']), 'fig')
    %end
    
    modelIndVarsMat{indL} = modelIndVars_mat;
    modelCoefMat{indL} = modelCoef_mat;
    Jw0s{indL} = Jw0;
    Jw1s{indL} = Jw1;
    
    %Pw1s{indL} = Pw1;
    %Qw0s{indL} = Qw0;
    %resiMap{indL} = resi_mat;
    %yhatFMap{indL} = yhatF_mat;
    %yhatRMap{indL} = yhatR_mat;
    for direction = 1:4  % order is LRUD.
        winFits{direction}{indL} = winFits_mat{direction};
    end
    
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    close(fmdlSlctn{indL})
end

%%   modelCoefMat

figwinCoefMat = cell(layerMax, 1);
%xtickvec = [(1:2*(1+2*wlagMax+2)) * tlagMax, 2*(1+2*wlagMax+2)*tlagMax + (1:(1+2*wlagMaxReg+2)) * tlagMaxReg];
xtickvec = [(1:5)*tlagMax];

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
    %avgArmaxLag = round(mean(winch2Pvec{indL}, 'omitnan'), 0);
    %avgRegLag = round(mean(winch1Pvec{indL}, 'omitnan'), 0);
    
    avgMaxOrder = [Jw0s{indL}, Jw1s{indL}];
    
    hold on
    for i=xtickvec
        line([i+0.5,i+0.5], [1,wmax])
    end
    
    title({['AR/adjWindow features chosen by IC - '];[ fsaveName0, '-', num2str(indL), 'L']; ...
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
save(fullfile(figuresDir, [fsaveName0, '_chosenOrder_Layers.mat']), 'Jw0s', 'Jw1s')
save(fullfile(figuresDir, [fsaveName0, '_winFits_drctn_Layers.mat']), 'winFits')


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
    %coefAR(xtickvec(5)+1:xtickvec(end)) = NaN;
    bar(coefAR, 'k');
    hold on
   
    %    bar(coefCurve)
    h = refline([0,0]); h.Color = [.3 .3 .3];
    set(gca, 'XTick', xtickvec);
    for i = xtickvec
        l = line([i+0.5,i+0.5], ylim); l.Color = [.3 .3 .3];
    end
    
    %avgArmaxLag = round(mean(winch2Pvec{indL}, 'omitnan'), 0);
    %avgRegLag = round(mean(winch1Pvec{indL}, 'omitnan'), 0);
    avgMaxOrder = [Jw0s{indL}, Jw1s{indL}];
    
    title({['AR/adjWindow features chosen by IC - '];[ fsaveName0, '-', num2str(indL), 'L']; ...
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

%%


%% plot pcorriCurve

%%  plot orders by AIC: AR, chan2
%%  plot orders by AIC: Reg, chan1

%% Chosen fit, Ftest, FPval, Rsq, cnum, sigma, resiMap
% pval, partialCorr, resiMap, Rsq, cnum, normality, ~, ~
%% chosen fit
%% saveas winFits

save(fullfile(figuresDir, ['GC_', fsaveName0, '_winFits_Layers.mat']), 'winFits')


%% draw outcome
%% draw ECDF of Fpval, LRTpval, normalityTest

figecdf = cell(layerMax, 4);
directionName = {'Left', 'Right', 'Outside', 'Inside'};
FpvecLyrDrctn = cell(layerMax, 4);
FDRLyrDrctn = cell(layerMax, 4);

% FpArr: 1l-LROI, 2l-LROI, ...
FpArr = nan(wmax * 4, layerMax);
for indL = 1:layerMax
    for direction = 1:4     % LRUD
        Fpvec = winFits{direction}{indL}(:, 8);
        FpvecLyrDrctn{indL, direction} = Fpvec;
        FpArr(1 + (direction - 1)*wmax:direction*wmax, indL) = Fpvec;
    end
end
FpAll = reshape(FpArr, [], 1);
fdrAll = mafdr(FpAll, 'BHFDR', true);       % mafdr adjustment
fdrArr = reshape(fdrAll, [], layerMax);

for indL = 1:layerMax
    for direction = 1:4     % LRUD
        figecdf{indL, direction} = figure('Visible', p.figFlag);
        
        Fpvec = winFits{direction}{indL}(:, 8);
        %LRTpvec = winFits{direction}{indL}(:, 10);
        KSpvec = winFits{direction}{indL}(:, 5);
        LBpvec = winFits{direction}{indL}(:, 12);
        % 202/01/10
        fdrbh = fdrArr(1 + (direction - 1)*wmax:direction*wmax, indL);
        FpvecLyrDrctn{indL, direction} = Fpvec;
        FDRLyrDrctn{indL, direction} = fdrbh;
        
        if ~all(isnan(Fpvec))
            [f,x] = ecdf(Fpvec);
            stairs(x, f)%, 'LineWidth', 2)
          
            hold on
            
            %[f,x] = ecdf(LRTpvec);
            %stairs(x, f, ':')%, 'LineWidth', 2)
            [f,x] = ecdf(KSpvec);
            stairs(x, f)
            [f,x] = ecdf(LBpvec);
            stairs(x, f)
            [f,x] = ecdf(fdrbh);
            stairs(x, f)
        else
            x = nan; f=nan; stairs(x,f)
        end
        
        legend('F-test', 'resiNormality', 'LB-test', 'BH-FDR', 'Location','southeast','Orientation','vertical');
        
        h=refline([1,0]); %h.Color='k';
        h=line([0.05, 0.05], [0,1]); h.Color='red';
        %h=refline([0, 0.5]);h.Color=[0.5 0.5 0.5];
        h=line([0.2, 0.2], [0,1]); h.Color='m';
        
        winProp = sum(Fpvec < 0.05)/sum(~isnan(Fpvec));
        medFPval = median(Fpvec, 'omitnan');
        
        subtitle0 = ['medianPval: ', num2str(round(medFPval, 3)), ', prop of <5%: ', num2str(round(winProp,2))];
        title0 = ['GC P-values, ', chan2Name, ' from ', directionName{direction}, ...
            ' to Center', ' -', num2str(indL), 'L'];
        title({title0; subtitle0})
        
        xlabel('P-value')
        ylabel('Cumulative Relative Frequency')
        
        ax = gca; ax.FontSize = 13;
        
        % save
        saveas(figecdf{indL, direction}, fullfile(figuresDir, ...
            [fsaveName0,'_from', directionName{direction}, '_', num2str(indL),'L_GCecdf.png']), 'png')
    end
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    for direction = 1:4
        close(figecdf{indL, direction})
    end
end

%% save
save(fullfile(figuresDir, [fsaveName0, '_Fpvec_FDR.mat']), ...
    'FpvecLyrDrctn', 'FDRLyrDrctn', 'FpArr', 'fdrArr')

%%  partial Rsquare

hist_rsq = cell(layerMax, direction);

for indL = 1:layerMax
    for direction = 1:4
        
        Fpvec = winFits{direction}{indL}(:, 8);
        partR2 = 100 * winFits{direction}{indL}(:, 11);
        ind0 = (Fpvec >= 0.05);
        ymaxR2 = max(10, max(partR2));
        %
        partR2sig = partR2;
        partR2sig(ind0) = NaN;
        %
        fdrbh = fdrArr(1 + (direction - 1)*wmax:direction*wmax, indL);
        indFDR = (fdrbh >= 0.2);
        partR2sigFDR = partR2;
        partR2sigFDR(indFDR) = NaN;
        
        hist_rsq{indL, direction} = figure('Visible', p.figFlag);
        %histogram(partR2, 'BinMethod', 'fd')
        plot(1:wmax, partR2)
        ylim([0, ymaxR2])
        
        hold on
        scatter(1:wmax, partR2sig, [], 'r')
        scatter(1:wmax, partR2sigFDR, [], '+')
        refline([0, mean(partR2, 'omitnan')]);
        title([fsaveName0, ' from ', directionName{direction}, ...
            ' to Center', ' -', num2str(indL), 'L', ' mean: ', num2str(mean(partR2, 'omitnan'))])
        xlabel('Window')
        ylabel('Partial R-square (%)')
        
        
        ax = gca; ax.FontSize = 10;
        
        %%
        saveas(hist_rsq{indL, direction}, fullfile(figuresDir, ...
            [fsaveName0,'_from', directionName{direction}, '_', num2str(indL),'L_partialRsquares.png']), 'png')
    end
end

%% close figures
pause(0.5)
for indL = 1:layerMax
    for direction = 1:4
        close(hist_rsq{indL, direction})
    end
end

%%  condition number

hist_cnum = cell(1, layerMax);

for indL = 1:layerMax
    
    %Fpvec = winFits{indL}(:, 8);
    cnumVec = winFits{1}{indL}(:, 2);
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
disp('====End of MD_GC_SPAR1ch_informationFlow ====')

end

