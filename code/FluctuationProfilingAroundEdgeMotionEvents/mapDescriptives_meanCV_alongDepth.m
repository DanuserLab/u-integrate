function mapDescriptives_meanCV_alongDepth(MD, iChan, maxLayer, chanName, ...
        figuresDir, varargin)
% mapDescriptives_meanCV_alongDepth Compute average or coeff. variation of 
% activities over time in each window, and Visualize them spatially and
% along the distance from edge.
%
% Updated:
% Jungsik Noh, 2018/01/30.
% Update by Qiongjing (Jenny) Zou, Aug 2018 

ip = inputParser;
ip.addParameter('figFlag', 'off');
ip.addParameter('impute', false);  
ip.addParameter('WithN', false);
ip.addParameter('omittedWindows', []);
ip.addParameter('Folding', false);
ip.addParameter('subFrames', []);
ip.addParameter('derivative', false);
ip.addParameter('topograph', 'on');
ip.addParameter('movingAvgSmoothing', true); % here, goal is clustering

ip.parse(varargin{:});
p = ip.Results;

%figFlag = p.figFlag;


%%  figuresDir setup
% figuresDir = fullfile(outDir, figDirName)           %% input
if ~isdir(figuresDir); mkdir(figuresDir); end

tmptext = ['mapDescriptives_meanSD_alongDepth_', 'inputParser.mat'];
save(fullfile(figuresDir, tmptext), 'p')


%%  getting Maps from channels

disp(chanName)


[fname0, MDpixelSize_, MDtimeInterval_, wmax, tmax, rawActmap, actmap_outl, imActmap] ...
        = mapOutlierImputation(MD, iChan, maxLayer, 'impute', p.impute, ...
            'WithN', p.WithN, 'omittedWindows', p.omittedWindows, ...
            'Folding', p.Folding, 'subFrames', p.subFrames, 'movingAvgSmoothing', p.movingAvgSmoothing); 
        

        
disp(['== MDpixelSize_: ', num2str(MDpixelSize_), ' =='])
disp(['== MDtimeInterval_: ', num2str(MDtimeInterval_), ' =='])


%%  .txt (export comma delimited files)

for indL = 1:maxLayer
    dlmwrite(fullfile(figuresDir, [fname0, '_', num2str(indL), 'L_actmap_outl.txt']), ...
                    actmap_outl{indL}, 'precision', 8)
    dlmwrite(fullfile(figuresDir, [fname0, '_', num2str(indL), 'L_imActmap.txt']), ...
                    imActmap{indL}, 'precision', 8)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Plots 


%% topomap topographMD: mean

meanlayers = cell(1, maxLayer);

if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    meanlayers{indL} = mean(actmap_outl{indL}, 2, 'omitnan');
    topoMap(:, indL) = meanlayers{indL};
end

title0 = [chanName, '-mean'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topograph_mean_', fname0, '.png']), 'png')
saveas(topomapFig, fullfile(figuresDir, ['/topograph_mean_', fname0, '.fig']), 'fig')

end

%%
save(fullfile(figuresDir, [fname0, '_', 'meanlayers.mat']), 'meanlayers')


%% histogram

histMeans = figure('Visible', p.figFlag);
tmp = cell2mat(meanlayers);
histogram(tmp(:), 'BinMethod', 'fd');
ylabel('Frequency')
xlabel('Intensity')
title(title0)
    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/hist_mean_', fname0, '.png']), 'png')

%%  BoxPlot per layer

bpmat = cell2mat(meanlayers);
BPlayer = figure('Visible', p.figFlag); 
boxplot(bpmat);
title(title0)
xlabel('Layers')
    ax = gca;
    ax.FontSize = 15;

saveas(BPlayer, fullfile(figuresDir, [fname0, '_mean_BPlayer.png']), 'png')

%%

%% topomap topographMD: SD

SDlayers = cell(1, maxLayer);

if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    SDlayers{indL} = std(actmap_outl{indL}, [], 2, 'omitnan');
    topoMap(:, indL) = SDlayers{indL};
end

title0 = [chanName, '-SD'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topograph_SD_', fname0, '.png']), 'png')
saveas(topomapFig, fullfile(figuresDir, ['/topograph_SD_', fname0, '.fig']), 'fig')

end

%%
save(fullfile(figuresDir, [fname0, '_', 'SDlayers.mat']), 'SDlayers')


%% histogram

histMeans = figure('Visible', p.figFlag);
tmp = cell2mat(SDlayers);
histogram(tmp(:), 'BinMethod', 'fd');
ylabel('Frequency')

title(title0)
    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/hist_SD_', fname0, '.png']), 'png')

%%  BoxPlot per layer

bpmat = cell2mat(SDlayers);
BPlayer = figure('Visible', p.figFlag); 
boxplot(bpmat);
title(title0)
xlabel('Layers')
    ax = gca;
    ax.FontSize = 15;

saveas(BPlayer, fullfile(figuresDir, [fname0, '_SD_BPlayer.png']), 'png')



%%
%% topomap topographMD: CV

CVlayers = cell(1, maxLayer);

if strcmp(p.topograph, 'on')

iWinProc = MD.getProcessIndex('WindowingProcess',1,0);
nBandMax_ = MD.processes_{iWinProc}.nBandMax_;
%
topoMap = nan(wmax, nBandMax_);

for indL = 1:maxLayer
    CVlayers{indL} = SDlayers{indL} ./ meanlayers{indL};
    topoMap(:, indL) = CVlayers{indL};
end

title0 = [chanName, '-CV'];

topomapFig = topographMD(MD, tmax, 1, topoMap, title0, p.figFlag);

ax = gca;
    ax.FontSize = 15;

%%
saveas(topomapFig, fullfile(figuresDir, ['/topograph_CV_', fname0, '.png']), 'png')
saveas(topomapFig, fullfile(figuresDir, ['/topograph_CV_', fname0, '.fig']), 'fig')

end

%%
save(fullfile(figuresDir, [fname0, '_', 'CVlayers.mat']), 'CVlayers')


%% histogram

histMeans = figure('Visible', p.figFlag);
tmp = cell2mat(CVlayers);
histogram(tmp(:), 'BinMethod', 'fd');
ylabel('Frequency')

title(title0)
    ax = gca;
    ax.FontSize = 15;

saveas(histMeans, fullfile(figuresDir, ['/hist_CV_', fname0, '.png']), 'png')

%%  BoxPlot per layer

bpmat = cell2mat(CVlayers);
BPlayer = figure('Visible', p.figFlag); 
boxplot(bpmat);
title(title0)
xlabel('Layers')
    ax = gca;
    ax.FontSize = 15;

saveas(BPlayer, fullfile(figuresDir, [fname0, '_CV_BPlayer.png']), 'png')

%%

%%
disp('====End of mapDescriptives_meanCV_alongDepth ====')
    
end


