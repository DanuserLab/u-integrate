function MLsummary_alongDepth(ML, iChan, maxLayer, chanName, ...
        analFolderName, varargin)
% MLsummary_alongDepth Summarize the output from
% mapDescriptives_meanCV_alongDepth().
%
% J Noh, 2018/01/30.



% load ML, example MD
ML.getMovies()
md1 = ML.getMovie(1);
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

%%
MDs = ML.getMovies();
num = numel(MDs);

chfname = ['Chan', num2str(iChan)];

ip = inputParser; 
ip.addParameter('outDirName', ['topolayers_', chfname]);
ip.addParameter('timeInterval', md1.timeInterval_);

parse(ip, varargin{:})
p = ip.Results;

fname_CV= [chfname, '_CVlayers.mat'];
fname_mean= [chfname, '_meanlayers.mat'];

fname_CVtopo= ['topograph_CV_', chfname, '.png'];
fname_meantopo= ['topograph_mean_', chfname, '.png'];

MDtimeIntvl = p.timeInterval;

iWinProc = md1.getProcessIndex('WindowingProcess',1,0);
winPerpSize = md1.processes_{iWinProc}.funParams_.PerpSize;
depthUM = md1.pixelSize_ * winPerpSize / 10^3;


%analNameAcf = 'mapDescriptives_0516_Imp0L4';  
%analNameXcf = 'mapCrossCorr_0510_Imp0L3';  
%% setting up parameters

outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end


%%  CVlayers

MDs = ML.getMovies();
num = numel(MDs);
%num = num    %%%  input
cellLabels = cell(num, 1);

CVArr = nan(num, maxLayer);

%
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    cellName = [cellLab0(1:end)];    % Generic cellName

    %%%% input  CVlayers (cell)
    load(fullfile(mdDir, analFolderName, fname_CV));   

    %%%%
    
    cellLabels{i} = cellName;
    copyfile(fullfile(mdDir, analFolderName, fname_CVtopo), ...
        fullfile(outDir, [cellLabels{i}, '_', fname_CVtopo]) )    
    
    %xcorrMat = Avg_autocor;
        %xcmean = mean(xcorrMat_tmp{indL}, 1, 'omitnan');
        %xcmean = Avg_autocor;
    
    tmp = cell2mat(CVlayers(1:maxLayer));
    CVvec = nanmean(tmp, 1);
        %xcmean = Avg_autocorLayers{1};
        %tmplen = numel(xcmean);
        %if (numel(xcmean) < lagSizeAcf)
        %    xcmeanMiddle = [xcmean, nan(1, lagSizeAcf - numel(xcmean))];
        %else
        %    xcmeanMiddle = xcmean(1:lagSizeAcf);
        %end
            
    CVArr(i, :) = CVvec;
    
end
  

%%
xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fCVlayers = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(CVArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(CVArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', CVArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5c .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    title([chanName, '-CV'])
    xlabel('Distance (um)');ylabel('Coef. of Variation')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];
    ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];

%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  


    

%% saveas
saveas(fCVlayers, fullfile(outDir, ['fCVlayers', '.png']), 'png')
save(fullfile(outDir, 'CVArr.mat'), 'CVArr')

%%
xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fCVlayers2 = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(CVArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(CVArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', CVArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5c .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    %title([chanName, '-CV'])
    %xlabel('Distance (um)');ylabel('Coef. of Variation')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    %ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];
    ax.XLim = [0, xx2(end)*1.1];
    %figure(fCVlayers);
    %title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)

saveas3format(fCVlayers2, outDir, 'fCVlayers2')


%%  meanlayers

MDs = ML.getMovies();
num = numel(MDs);
%num = num    %%%  input
cellLabels = cell(num, 1);

meanArr = nan(num, maxLayer);
normMeanArr = nan(num, maxLayer);
baseMeanArr = nan(num, maxLayer);

%
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    cellName = [cellLab0(1:end)];    % Generic cellName

    %%%% input  CVlayers (cell)
    load(fullfile(mdDir, analFolderName, fname_mean));   

    %%%%
    
    cellLabels{i} = cellName;
    copyfile(fullfile(mdDir, analFolderName, fname_meantopo), ...
        fullfile(outDir, [cellLabels{i}, '_', fname_meantopo]) )    
    
    %xcorrMat = Avg_autocor;
        %xcmean = mean(xcorrMat_tmp{indL}, 1, 'omitnan');
        %xcmean = Avg_autocor;
    
    tmp = cell2mat(meanlayers(1:maxLayer));
    meanvec = nanmean(tmp, 1);
        %xcmean = Avg_autocorLayers{1};
        %tmplen = numel(xcmean);
        %if (numel(xcmean) < lagSizeAcf)
        %    xcmeanMiddle = [xcmean, nan(1, lagSizeAcf - numel(xcmean))];
        %else
        %    xcmeanMiddle = xcmean(1:lagSizeAcf);
        %end
            
    meanArr(i, :) = meanvec;
    normMeanArr(i, :) = meanvec ./ meanvec(1);
    baseMeanArr(i, :) = meanvec ./ meanvec(end);
    
end
  

%%
xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fmeanlayers = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(normMeanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(normMeanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', normMeanArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    title([chanName, '-mean'])
    xlabel('Distance (um)');ylabel('Normalized activity')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];    
    ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];

%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  


    %fCVlayers2 = figure();    %'Visible', p.figFlag);   %%  name


%% saveas
saveas(fmeanlayers, fullfile(outDir, ['fmeanlayers', '.png']), 'png')
save(fullfile(outDir, 'normMeanArr.mat'), 'normMeanArr')

%%
xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fmeanlayers2 = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(normMeanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(normMeanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', normMeanArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    %title([chanName, '-mean'])
    %xlabel('Distance (um)');ylabel('Normalized activity')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    %ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];
    ax.XLim = [0, xx2(end)*1.1];
%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    %figure(fmeanlayers);
    
    %title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)

saveas3format(fmeanlayers2, outDir, 'fmeanlayers2')


%%
%%  meanlayers0

xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fmeanlayers0 = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(meanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(meanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', meanArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    title([chanName, '-mean'])
    xlabel('Distance (um)');ylabel('Averaged activity')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];    
    ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];

%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  


    %fCVlayers2 = figure();    %'Visible', p.figFlag);   %%  name


%% saveas
saveas(fmeanlayers0, fullfile(outDir, ['fmeanlayers0', '.png']), 'png')
save(fullfile(outDir, 'meanArr.mat'), 'meanArr')

%%
xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fmeanlayers02 = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(meanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(meanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', meanArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    %title([chanName, '-mean'])
    %xlabel('Distance (um)');ylabel('Normalized activity')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    %ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];
    ax.XLim = [0, xx2(end)*1.1];
%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    %figure(fmeanlayers);
    
    %title(''); xlabel(''); ylabel('')
    ax = gca; ax.FontSize = 25;
    pause(0.1)

saveas3format(fmeanlayers02, outDir, 'fmeanlayers02')


%%  2018/05/03
%%  meanlayer_baseline

xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    fmeanlayersBase = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(baseMeanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(baseMeanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    plot(xx2', baseMeanArr', 'Color', 'k')
    %h1 = refline([0,0]); h1.Color = [.5 .5 .5];
    %h = line([lagMax+1, lagMax+1], ylim); h.Color = [.5 .5 .5];
    %hold on
    %p1 = plot(xx', nanmean(CVArr, 1)', 'ro-', 'MarkerSize', 10);
    %p1.LineWidth = 2;
    
    title([chanName, '-mean'])
    xlabel('Distance (um)');ylabel('Normalized activity')
    set(gca, 'XGrid', 'on')
    
    h=refline(0, 1);h.Color = [.5 .5 .5];
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];    
    ax.FontSize = 14;
    %ax.XLim = [0, maxLayer+1];

%%%%  !!!!!!!!!!!!!!!!!!!!!!!!!!!!  


    %fCVlayers2 = figure();    %'Visible', p.figFlag);   %%  name


%% saveas
saveas(fmeanlayersBase, fullfile(outDir, ['fmeanlayersBase', '.png']), 'png')
save(fullfile(outDir, 'baseMeanArr.mat'), 'baseMeanArr')

%%  paper plot

xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    ffmeanlayersBase = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(baseMeanArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(baseMeanArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
   
    title([chanName, '-mean'])
    xlabel('Distance (um)');ylabel('Normalized activity')
    set(gca, 'XGrid', 'on')
    
    h=refline(0, 1);h.Color = 'k';
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];    
    ax.FontSize = 14;

saveas(ffmeanlayersBase, fullfile(outDir, ['ffmeanlayersBase', '.png']), 'png')
saveas(ffmeanlayersBase, fullfile(outDir, ['ffmeanlayersBase', '.fig']), 'fig')

%%
figure(ffmeanlayersBase)
title(''); xlabel(''); ylabel('');
    ax = gca; ax.FontSize = 25;
    pause(0.1)

saveas3format(ffmeanlayersBase, outDir, 'ffmeanlayersBase2')

%%  paper plot2

xx = 1:maxLayer; 
xx2 = xx .* depthUM;

    ffCVlayers = figure();    %'Visible', p.figFlag);   %%  name
    
    sem = std(CVArr, [], 1, 'omitnan')./sqrt(num);
    
    s1 = shadedErrorBarV2(xx2, nanmean(CVArr, 1), 2*sem, 'lineprops', ...
        {'ro-', 'MarkerSize', 10});
    s1.mainLine.LineWidth = 2;    
    
    hold on
    
    title([chanName, '-CV'])
    xlabel('Distance (um)');ylabel('Coef. of Variation')
    set(gca, 'XGrid', 'on')
    
    ax = gca;
    ax.XLim = [0, xx2(end)*1.1];
    ax.FontSize = 14;

saveas(ffCVlayers, fullfile(outDir, ['ffCVlayers', '.png']), 'png')
saveas(ffCVlayers, fullfile(outDir, ['ffCVlayers', '.fig']), 'fig')

%%
figure(ffCVlayers)
title(''); xlabel(''); ylabel('');
    ax = gca; ax.FontSize = 25;
    pause(0.1)

saveas3format(ffCVlayers, outDir, 'ffCVlayers2')




%%
disp('==== MLsummary_alongDepth is finished!! ====')
disp('== :) close all ')

end
