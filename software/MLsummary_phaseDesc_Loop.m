function MLsummary_phaseDesc_Loop(md1, ch0ActmapName, maxLayer, outDir, ...
                                analName, MDs, num, ch1ActmapName, frSize, MDtimeInterval, figFlag)
% MLsummary_phaseDesc_Loop
%
% Updated:
% J Noh, 2018/10/29. Add a 'figFlag' option.
% J Noh, 2018/02/22. Add a 'MDtimeInterval' argument to handle the 'Folding'
% option.
% J Noh, 2018/01/30
%
% Copyright (C) 2022, Danuser Lab - UTSouthwestern 
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
                            
                            
%MDs = ML.getMovies();
%num = numel(MDs);

%
num = num    %%%  input

cellLabels = cell(num, 1);

MLch0Zsamples = cell(num, 1);
MLch1Zsamples = cell(num, 1);
 
 
for i = 1:num
    md = MDs{i};
    mdDir = md.outputDirectory_;
    [folderName, cellLab0, ~] = fileparts(mdDir);
    %cellName = [folderName((end-4):end), '_', cellLab0];   %% Adjust name extraction
    cellName = [cellLab0(1:end)];   %% Adjust name extraction
    %[~, cellName] = fileparts(folderName);

    %%%% input
    load(fullfile(mdDir, analName, [ch1ActmapName, '-ch0Zsamples.mat']));   % xcorrMat_tmp
    load(fullfile(mdDir, analName, [ch1ActmapName, '-ch1Zsamples.mat']));   % permXcorrMean
    %%%%
    
    cellLabels{i} = cellName;
    for indL = 1:maxLayer
        fnametmp = [ch1ActmapName, '_meanTS_', num2str(indL), 'L.png'];
        copyfile(fullfile(mdDir, analName, fnametmp), ...
            fullfile(outDir, [cellLabels{i}, '_', fnametmp]) )    
    end
    
    MLch0Zsamples{i} = ch0Zsamples;
    MLch1Zsamples{i} = ch1Zsamples;
end

%
save(fullfile(outDir, [ch1ActmapName, '-MLch01Zsamples.mat']), ...
      'cellLabels', 'MLch0Zsamples', 'MLch1Zsamples', 'ch1ActmapName', ...
      'ch0ActmapName', 'maxLayer', 'analName')


%%  Integrate information using only mean vectors per movie
% 2017/12/11, ch0Zsamples, frSize

headSize = (size(ch0Zsamples{1}{1}, 2) - frSize) / 2;

samplingFrSize = frSize; 
samplingBw = (samplingFrSize-1)/2;

% Onset analysis: plot sampled Zscores 

timeAxis = (-samplingBw:1:samplingBw)*MDtimeInterval;

for indL = 1:maxLayer
    
    mat1M = [];
    for k = 1:num
        tmp = cell2mat(MLch0Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');
        mat1M = [mat1M; tmpM];
    end
    
    mat2M = [];
    for k = 1:num
        tmp = cell2mat(MLch1Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');        
        mat2M = [mat2M; tmpM];
    end
    
    %mat1tmp = mat1(~any(isnan(mat1), 2), :);
    %mat2tmp = mat2M(~any(isnan(mat2M), 2), :);
    sem1 = std(mat1M, [], 1)./sqrt(num);
    sem2 = std(mat2M, [], 1)./sqrt(num);    

    meanTS_fig{indL} = figure('Visible', figFlag); 
    s1 = shadedErrorBarV2(timeAxis, mean(mat1M, 1, 'omitnan'), 2*sem1, 'lineprops', '-r');
    %s1 = plot(timeAxis, mean(mat1M, 1, 'omitnan'), '-r', 'LineWidth', 2);
    hold on
    s2 = shadedErrorBarV2(timeAxis, mean(mat2M, 1, 'omitnan'), 2*sem2, 'lineprops', '-b');
    %s2 = plot(timeAxis, mean(mat2M, 1, 'omitnan'), '-b', 'LineWidth', 2);
    s1.mainLine.LineWidth = 2;
    s2.mainLine.LineWidth = 2;
    
    %p1 = plot(timeAxis, mat1M, '-r', 'LineWidth', 0.1);
    %p2 = plot(timeAxis, mat2M, '-b', 'LineWidth', 0.1);
    
    title(['Average of locally sampled standardized TS -', num2str(indL), 'L']);
    %title0 = [num2str(indL), 'L, bandwidth (lag) around reference timing: ', num2str(samplingBw)];
    %title({title1; title0})
    
    xlabel('Time (s)');ylabel('Standardized TS')
    
    hold on; 
    refline([0, 0])
    h = line([0 0], ylim);
    
    legend([s1.mainLine, s2.mainLine], {ch0ActmapName, ch1ActmapName}, 'Location', 'northwest')    
end



%% saveas
    for indL = 1:maxLayer
        saveas(meanTS_fig{indL}, fullfile(outDir, [ch1ActmapName, ...
            '_MLZsamples_meanTS_',  num2str(indL), 'L.png']), 'png')
    end


    for indL = 1:maxLayer
        saveas(meanTS_fig{indL}, fullfile(outDir, [ch1ActmapName, ...
            '_MLZsamples_meanTS_',  num2str(indL), 'L.fig']), 'fig')
    end    



%%  Visualize information using only mean vectors per movie

headSize = (size(ch0Zsamples{1}{1}, 2) - frSize) / 2;

samplingFrSize = frSize; 
samplingBw = (samplingFrSize-1)/2;

% Onset analysis: plot sampled Zscores 

timeAxis = (-samplingBw:1:samplingBw)*md1.timeInterval_;

for indL = 1:maxLayer
    
    mat1M = [];
    for k = 1:num
        tmp = cell2mat(MLch0Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');
        mat1M = [mat1M; tmpM];
    end
    
    mat2M = [];
    for k = 1:num
        tmp = cell2mat(MLch1Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');        
        mat2M = [mat2M; tmpM];
    end
    
    %mat1tmp = mat1(~any(isnan(mat1), 2), :);
    %mat2tmp = mat2M(~any(isnan(mat2M), 2), :);
    %sem1 = std(mat1M, [], 1)./sqrt(num);
    %sem2 = std(mat2M, [], 1)./sqrt(num);    

    summaryTS_fig{indL} = figure('Visible', figFlag); 
    
    s1 = plot(timeAxis, mean(mat1M, 1, 'omitnan'), '-r', 'LineWidth', 2);
    hold on
    
    s2 = plot(timeAxis, mean(mat2M, 1, 'omitnan'), '-b', 'LineWidth', 2);
    
    p1 = plot(timeAxis, mat1M, '-r', 'LineWidth', 0.5);
    p2 = plot(timeAxis, mat2M, '-b', 'LineWidth', 0.5);
    
    title(['Average of locally sampled standardized TS -', num2str(indL), 'L']);
    %title0 = [num2str(indL), 'L, bandwidth (lag) around reference timing: ', num2str(samplingBw)];
    %title({title1; title0})
    
    xlabel('Time (s)');ylabel('Standardized TS')
    
    refline([0, 0])
    hold on; h = line([0 0], ylim);
    
    
    legend([s1, s2], {'Vel', ch1ActmapName}, 'Location', 'northwest')    
end


%% saveas
    for indL = 1:maxLayer
        saveas(summaryTS_fig{indL}, fullfile(outDir, [ch1ActmapName, ...
            '_MLZsamples_summaryTS_',  num2str(indL), 'L.png']), 'png')
    end


    for indL = 1:maxLayer
        saveas(summaryTS_fig{indL}, fullfile(outDir, [ch1ActmapName, ...
            '_MLZsamples_summaryTS_',  num2str(indL), 'L.fig']), 'fig')
    end    

    

%%  means fig2

for indL = 1:maxLayer
    
    mat1M = [];
    for k = 1:num
        tmp = cell2mat(MLch0Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');
        mat1M = [mat1M; tmpM];
    end
    
    mat2M = [];
    for k = 1:num
        tmp = cell2mat(MLch1Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');        
        mat2M = [mat2M; tmpM];
    end
    
    %mat1tmp = mat1(~any(isnan(mat1), 2), :);
    %mat2tmp = mat2M(~any(isnan(mat2M), 2), :);
    sem1 = std(mat1M, [], 1)./sqrt(num);
    sem2 = std(mat2M, [], 1)./sqrt(num);    

    meanTS_fig2{indL} = figure('Visible', figFlag); 
    s1 = shadedErrorBarV2(timeAxis, mean(mat1M, 1, 'omitnan'), 2*sem1, 'lineprops', '-r');
    %s1 = plot(timeAxis, mean(mat1M, 1, 'omitnan'), '-r', 'LineWidth', 2);
    hold on
    s2 = shadedErrorBarV2(timeAxis, mean(mat2M, 1, 'omitnan'), 2*sem2, 'lineprops', '-b');
    %s2 = plot(timeAxis, mean(mat2M, 1, 'omitnan'), '-b', 'LineWidth', 2);
    s1.mainLine.LineWidth = 2;
    s2.mainLine.LineWidth = 2;
    
    ax = gca; ax.FontSize = 25;
    
    hold on; 
    refline([0, 0])
    h = line([0 0], ylim); 
    
        pause(0.1)
        
        saveas3format(meanTS_fig2{indL}, outDir, ['paper1_', ch1ActmapName, '_MLZsamples_meanTS_',  num2str(indL), 'L'])
        pause(0.1)
%    end    

end    


%%  means fig3

for indL = 1:maxLayer
    
    mat1M = [];
    for k = 1:num
        tmp = cell2mat(MLch0Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');
        mat1M = [mat1M; tmpM];
    end
    
    mat2M = [];
    for k = 1:num
        tmp = cell2mat(MLch1Zsamples{k}{indL});
        % 2017/12/11, ch0Zsamples, frSize
        tmp = tmp(:, headSize+1:end-headSize);
        tmpM = mean(tmp, 1, 'omitnan');        
        mat2M = [mat2M; tmpM];
    end
    
    %mat1tmp = mat1(~any(isnan(mat1), 2), :);
    %mat2tmp = mat2M(~any(isnan(mat2M), 2), :);
    sem1 = std(mat1M, [], 1)./sqrt(num);
    sem2 = std(mat2M, [], 1)./sqrt(num);    

    meanTS_fig2{indL} = figure('Visible', figFlag); 
    s1 = shadedErrorBarV2(timeAxis, mean(mat1M, 1, 'omitnan'), 2*sem1, 'lineprops', '-r');
    %s1 = plot(timeAxis, mean(mat1M, 1, 'omitnan'), '-r', 'LineWidth', 2);
    hold on
    s2 = shadedErrorBarV2(timeAxis, mean(mat2M, 1, 'omitnan'), 2*sem2, 'lineprops', '-b');
    %s2 = plot(timeAxis, mean(mat2M, 1, 'omitnan'), '-b', 'LineWidth', 2);
    s1.mainLine.LineWidth = 2;
    s2.mainLine.LineWidth = 2;
    
    p2 = plot(timeAxis, mat2M, '-b', 'LineWidth', 0.5);    
    
    ax = gca; ax.FontSize = 25;
    
    hold on; 
    refline([0, 0])
    h = line([0 0], ylim); 
    
        pause(0.1)
        
        saveas3format(meanTS_fig2{indL}, outDir, ['paper2_', ch1ActmapName, '_MLZsamples_summary2_',  num2str(indL), 'L'])
        pause(0.1)
%    end    

end    

    
    

%%  EOF

end


