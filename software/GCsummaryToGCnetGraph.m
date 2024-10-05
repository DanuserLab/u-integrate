function GCsummaryToGCnetGraph(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayer)
% GCsummaryToGCnetGraph() Collect output from 
%
% Updates:
% J Noh, 2021/02/11. Rename and add comments.
% J Noh, 2020/01/04. Add GC-strengths in Network Diagrams using
% partial-R2.
% J Noh, 2019/08/06.
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


%% GC pval file names

outLab = cell2mat(chanDetailedNames);
chanDetailedNames = chanDetailedNames(:);

matfiles = dir(fullfile(GCsummaryDir, 'GC_Fpval*.mat'));
matfilenames = {matfiles.name};
matfilenames = matfilenames(:);
numfiles = numel(matfilenames);

%% first filter out irrelavant GC output files
indRelavant = false(numfiles, 1);

for k = 1:numfiles
    GCstrParts = strsplit(matfilenames{k}, {'_', '.mat'});
    stid = cellfun(@(x) find(strcmp(GCstrParts, x)), chanDetailedNames, 'UniformOutput', false);
    tmp = cellfun(@isempty, stid);
    indRelavant(k) = (sum(tmp == 0) >= 2);
end

matfilenames2 = matfilenames(indRelavant);
numfiles2 = numel(matfilenames2);


%% chanNames = {X;Y;Z}, match pval files to yx, zx, xy, zy, xz, yz
% directed edges are ordered by the order of elements in connectivity
% matrix, which is yx, zx, xy, zy, xz, yz.
% Then, rank of x,y,z in the file name will be 213,231,123,321,132,312.

GCidvec = nan(numfiles2, 1);

for k = 1:numfiles2
    % stid = cellfun(@(x) strfind(matfilenames{k}, x), chanDetailedNames, 'UniformOutput', false);
    GCstrParts = strsplit(matfilenames2{k}, {'_', '.mat'});
    stid = cellfun(@(x) find(strcmp(GCstrParts, x)), chanDetailedNames, 'UniformOutput', false);

    emptyid = cellfun(@isempty, stid);
    if any(emptyid)         %        stid2 = stid;
        stid{emptyid} = NaN;
    end
        stid3 = cell2mat(stid);
    [~, s] = sort(stid3); [~,r] = sort(s);
    %
    rInd = r(1)*100 + r(2)*10 + r(3);
    switch rInd
        case 213; GCind = 1;    % fr y to x given z
        case 231; GCind = 2;
        case 123; GCind = 3;
        case 321; GCind = 4;
        case 132; GCind = 5;
        case 312; GCind = 6;
    end
    GCidvec(k) = GCind;
end

matfilenamesSrtd = matfilenames2(GCidvec);
% mat indexing
cmrow = [2, 3, 1, 3, 1, 2];
cmcol = [1, 1, 2, 2, 3, 3];

%% computed one-sided rank test pvalues and link flags. make connectivity matrix

connMat = zeros(numel(chanDetailedNames), numel(chanDetailedNames), maxLayer);
partR2Mat = zeros(numel(chanDetailedNames), numel(chanDetailedNames), maxLayer);

for k = 1:numfiles2
    S = load(fullfile(GCsummaryDir, matfilenamesSrtd{k}));
    % mgdMedFpvalMat
    FpvalMat = S.medFpvalMat(:, 1:maxLayer);
    signRankPvec = nan(1, size(FpvalMat,2));
    for l = 1:size(FpvalMat, 2)
        pvals = FpvalMat(:, l);
        if any(~isnan(pvals))
            signRankPvec(l) = signrank(pvals, 0.05, 'tail', 'left');
        else
            signRankPvec(l) = NaN;
        end
    end
    signRankPvec1 = round(signRankPvec, 4);
    % later GC str: medmedFpvalVec = round(median(FpvalMat, 1, 'omitnan'), 4);
    % GC str by median partial R2
    medmedpartR2Mat = round(100*median(S.medpartR2Mat, 1, 'omitnan'), 0);
    
    sigIndic = (signRankPvec1 < 0.05)';
    for l = 1:size(FpvalMat, 2)
        connMat(cmrow(k), cmcol(k), l) = sigIndic(l);
        partR2Mat(cmrow(k), cmcol(k), l) = medmedpartR2Mat(l);
    end
end

connMatwt = connMat .* partR2Mat;

%% make directed graph, draw

grp = cell(maxLayer, 1);
fgrp = cell(maxLayer, 1);

for l = 1:maxLayer
    grp{l} = digraph(connMatwt(:,:,l), chanNamesforNet);
    fgrp{l} = figure;
    % Graph drawing by hierarchical graph drawing method using a
    % pre-assumption that edgeMotion is the end node.
    p = plot(grp{l}, 'layout', 'layered', ...
        'EdgeLabel', grp{l}.Edges.Weight, 'EdgeFontSize', 15, ...
        'Direction', 'down', 'Sinks', chanNamesforNet{3}, 'AssignLayers', 'alap');
    p.NodeColor = 'k'; p.Marker = 'o';
    p.MarkerSize = 8; p.NodeFontSize = 25;
    p.EdgeColor = 'r'; 
    % Edge scale: LW, 2~18~34 (0,50,100), Arrow, 20~60~100
    wtVec = grp{l}.Edges.Weight;
    arSz = 20 + (100-20) * wtVec ./ 100;
    lnWd = 2 + (34 - 2) * wtVec ./ 100;

    % to fix weired error from GUI run (2024)
    if ~isempty(arSz)
        if isfinite(arSz) 
            p.ArrowSize = arSz; p.LineWidth = lnWd;
        end
    end

    p.EdgeAlpha = 0.9;
    ax = gca;
    %ax.Position(3) = 0.6;
    ax.Position = [0.1300 0.1450 0.75 0.7850];
    
    %
    title0 = ['Granger-Causal Pathway - ', num2str(l), 'L'];
    title(title0)
end

%% save

outfname = ['GCnetworkGraph_Output_', outLab, '.mat'];
save(fullfile(GCsummaryDir, outfname), 'connMat', 'connMatwt', 'maxLayer', 'chanDetailedNames', 'chanNamesforNet')

pause(1)

for l = 1:maxLayer
    saveas(fgrp{l}, fullfile(GCsummaryDir, ['GCnetGraph_', num2str(l), 'L_', outLab, '.png']), 'png')
    %saveas(fgrp{l}, fullfile(GCsummaryDir, 'GCnetGraph_', num2str(l), 'L.fig'), 'fig')
    saveas3format(fgrp{l}, GCsummaryDir, ['GCnetGraph_', num2str(l), 'L_', outLab]);
end

%% make directed graph, draw2

grp = cell(maxLayer, 1);
fgrp = cell(maxLayer, 1);

for l = 1:maxLayer
    grp{l} = digraph(connMatwt(:,:,l), chanNamesforNet);
    fgrp{l} = figure;
    % Graph drawing by hierarchical graph drawing method using a
    % pre-assumption that edgeMotion is the end node.
    p = plot(grp{l}, 'layout', 'layered', ...
        'Direction', 'down', 'Sinks', chanNamesforNet{3}, 'AssignLayers', 'alap');
    p.NodeColor = 'k'; p.Marker = 'o';
    p.MarkerSize = 8; p.NodeFontSize = 25;
    p.EdgeColor = 'r'; 
    % Edge scale: LW, 2~18~34 (0,50,100), Arrow, 20~60~100
    wtVec = grp{l}.Edges.Weight;
    arSz = 20 + (100-20) * wtVec ./ 100;
    lnWd = 2 + (34 - 2) * wtVec ./ 100;

    % to fix weired error from GUI run (2024)
    if ~isempty(arSz)
        if isfinite(arSz) 
            p.ArrowSize = arSz; p.LineWidth = lnWd;
        end
    end
    
    p.EdgeAlpha = 0.9;
    ax = gca;
    %ax.Position(3) = 0.6;
    ax.Position = [0.1300 0.1450 0.75 0.7850];
    %
    title0 = ['Granger-Causal Pathway - ', num2str(l), 'L'];
    title(title0)
    %  %oldPos = fgrp{l}.PaperPosition(3:4);%newPos = 0.8 * oldPos; %fgrp{l}.PaperPosition(3:4) = newPos;
end

%% save2

pause(1)

for l = 1:maxLayer
    saveas(fgrp{l}, fullfile(GCsummaryDir, ['GCnetGraph2_', num2str(l), 'L_', outLab, '.png']), 'png')
    %saveas(fgrp{l}, fullfile(GCsummaryDir, 'GCnetGraph_', num2str(l), 'L.fig'), 'fig')
    saveas3format(fgrp{l}, GCsummaryDir, ['GCnetGraph2_', num2str(l), 'L_', outLab]);
end

%% make directed graph, draw3

grp = cell(maxLayer, 1);
fgrp = cell(maxLayer, 1);

for l = 1:maxLayer
    grp{l} = digraph(connMatwt(:,:,l), chanNamesforNet);
    fgrp{l} = figure;
    % Graph drawing by hierarchical graph drawing method using a
    % pre-assumption that edgeMotion is the end node.
    p = plot(grp{l}, 'layout', 'layered', ...
        'Direction', 'down', 'Sinks', chanNamesforNet{3}, 'AssignLayers', 'alap');
    p.NodeColor = 'k'; p.Marker = 'o';
    p.MarkerSize = 8; p.NodeFontSize = 25;
    p.EdgeColor = 'r'; 
    % Edge scale: LW, 2~18~34 (0,50,100), Arrow, 20~60~100
    wtVec = grp{l}.Edges.Weight;
    %arSz = 20 + (100-20) * wtVec ./ 100;
    %lnWd = 2 + (34 - 2) * wtVec ./ 100;
    p.ArrowSize = 25; p.LineWidth = 2;
    p.EdgeAlpha = 0.9;
    ax = gca;
    %ax.Position(3) = 0.6;
    ax.Position = [0.1300 0.1450 0.75 0.7850];
    %
    title0 = ['Granger-Causal Pathway - ', num2str(l), 'L'];
    title(title0)
    %  %oldPos = fgrp{l}.PaperPosition(3:4);%newPos = 0.8 * oldPos; %fgrp{l}.PaperPosition(3:4) = newPos;
end

%% save3

pause(1)

for l = 1:maxLayer
    saveas(fgrp{l}, fullfile(GCsummaryDir, ['GCnetGraph3_', num2str(l), 'L_', outLab, '.png']), 'png')
    %saveas(fgrp{l}, fullfile(GCsummaryDir, 'GCnetGraph_', num2str(l), 'L.fig'), 'fig')
    saveas3format(fgrp{l}, GCsummaryDir, ['GCnetGraph3_', num2str(l), 'L_', outLab]);
end

%%
disp('==== GCsummaryToGCnetGraph is done! ====') 

end