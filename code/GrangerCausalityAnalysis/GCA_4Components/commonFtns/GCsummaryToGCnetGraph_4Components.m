function GCsummaryToGCnetGraph_4Components(GCsummaryDir, chanDetailedNames, chanNamesforNet, maxLayer)
% GCsummaryToGCnetGraph_4Components() Collect GCA output and generate GC
% network diagrams.

%% === PRE-PROCESSING BLOCK (file parsing, GC direction, p-value thresholds) ===
outLab = cell2mat(chanDetailedNames);
chanDetailedNames = chanDetailedNames(:);
matfiles = dir(fullfile(GCsummaryDir, 'GC_Fpval*.mat'));
matfilenames = {matfiles.name}';
numfiles = numel(matfilenames);

indRelavant = false(numfiles, 1);
for k = 1:numfiles
    GCstrParts = strsplit(matfilenames{k}, {'_', '.mat'});
    stid = cellfun(@(x) find(strcmp(GCstrParts, x)), chanDetailedNames, 'UniformOutput', false);
    tmp = cellfun(@isempty, stid);
    indRelavant(k) = (sum(tmp == 0) >= 2);
end
matfilenames2 = matfilenames(indRelavant);
numfiles2 = numel(matfilenames2);

GCidvec = nan(numfiles2, 1);
availableIndices = 1:12;
equivalenceGroups = {[2134, 2143], [3124, 3142], [4123, 4132], [1234, 1243], ...
    [3214, 3241], [4213, 4231], [1324, 1342], [2314, 2341], ...
    [4312, 4321], [1423, 1432], [2413, 2431], [3412, 3421]};
for k = 1:numfiles2
    GCstrParts = strsplit(matfilenames2{k}, {'_', '.mat'});
    stid = cellfun(@(x) find(strcmp(GCstrParts, x)), chanDetailedNames, 'UniformOutput', false);
    emptyid = cellfun(@isempty, stid);
    stid(emptyid) = {NaN};
    stid3 = cell2mat(stid);
    [~, s] = sort(stid3); [~, r] = sort(s);
    rInd = r(1)*1000 + r(2)*100 + r(3)*10 + r(4);
    groupIndex = find(cellfun(@(group) any(ismember(group, rInd)), equivalenceGroups));
    GCidvec(k) = ~isempty(groupIndex) * availableIndices(1);
    if ~isempty(groupIndex), availableIndices(1) = []; end
end
matfilenamesSrtd = matfilenames2(GCidvec);

%%%%  mat indexing
cmrow =[1,1,1,2,2,2,3,3,3,4,4,4];
cmcol =[2,3,4,1,3,4,1,2,4,1,2,3];

%% Compute one-sided rank test p-values and make GC connectivity matrices with 5%/10% thresholds

connMat = zeros(numel(chanDetailedNames), numel(chanDetailedNames), maxLayer);
partR2Mat = zeros(numel(chanDetailedNames), numel(chanDetailedNames), maxLayer);
edgeColorMat = zeros(numel(chanDetailedNames), numel(chanDetailedNames), maxLayer);  % 0: none, 1: red, 2: gray

for k = 1:numfiles2
    S = load(fullfile(GCsummaryDir, matfilenamesSrtd{k}));
    FpvalMat = S.medFpvalMat(:, 1:maxLayer);
    signRankPvec = nan(1, size(FpvalMat, 2));
    
    for l = 1:size(FpvalMat, 2)
        pvals = FpvalMat(:, l);
        if any(~isnan(pvals))
            signRankPvec(l) = signrank(pvals, 0.05, 'tail', 'left');
        else
            signRankPvec(l) = NaN;
        end
    end
    
    signRankPvec1 = round(signRankPvec, 4);
    medR2 = round(100 * median(S.medpartR2Mat, 1, 'omitnan'), 0);
    
    for l = 1:maxLayer
        % Connectivity (strict 5% significance)
        if signRankPvec1(l) < 0.05
            connMat(cmrow(k), cmcol(k), l) = 1;
            edgeColorMat(cmrow(k), cmcol(k), l) = 1;  % Red edge
        elseif signRankPvec1(l) < 0.10
            connMat(cmrow(k), cmcol(k), l) = 0;       % Not significant at 5%
            edgeColorMat(cmrow(k), cmcol(k), l) = 2;  % Gray edge (moderate evidence)
        end
        
        partR2Mat(cmrow(k), cmcol(k), l) = medR2(l);  % Store R² regardless of significance
    end
end

connMatwt = connMat .* partR2Mat;  % Weighted matrix: strength only where p < 0.05


%% === Build plotMat: includes both red and gray edges ===
plotMat = zeros(size(connMat));
for i = 1:size(connMat,1)
    for j = 1:size(connMat,2)
        for l = 1:maxLayer
            if edgeColorMat(i,j,l) > 0  % red or gray
                plotMat(i,j,l) = partR2Mat(i,j,l);  % include weight for arrow size
            end
        end
    end
end
%% === PLOT STYLE 1: Threshold-scaled plot ===
grp = cell(maxLayer,1); fgrp = cell(maxLayer,1);
for l = 1:maxLayer
    connTemp = plotMat(:,:,l);  % use plotMat here
    colorTemp = edgeColorMat(:,:,l);
    grp{l} = digraph(connTemp, chanNamesforNet);
    fgrp{l} = figure;
    p = plot(grp{l}, 'layout','layered', 'Direction','down', ...
        'Sinks',chanNamesforNet{4}, 'AssignLayers','alap', ...
        'EdgeLabel', grp{l}.Edges.Weight, 'EdgeFontSize', 15);
    p.NodeColor = 'k'; p.Marker = 'o'; p.MarkerSize = 8; p.NodeFontSize = 25;
    wtVec = grp{l}.Edges.Weight;
    arSz = 20 + (100-20) * wtVec ./ 100;
    lnWd = 2 + (34 - 2) * wtVec ./ 100;
    p.ArrowSize = arSz; p.LineWidth = lnWd;
    p.EdgeAlpha = 0.9;

    edgeColors = zeros(numedges(grp{l}),3);
    for e = 1:numedges(grp{l})
        src = find(strcmp(chanNamesforNet, grp{l}.Edges.EndNodes{e,1}));
        tgt = find(strcmp(chanNamesforNet, grp{l}.Edges.EndNodes{e,2}));
        if colorTemp(src,tgt)==1
            edgeColors(e,:) = [1 0 0];         % red
        elseif colorTemp(src,tgt)==2
            edgeColors(e,:) = [0.5 0.5 0.5];   % gray
        end
    end
    p.EdgeColor = edgeColors;
    title(['Granger-Causal Pathway - ', num2str(l), 'L']);
    saveas(fgrp{l}, fullfile(GCsummaryDir, ['GCnetGraph_', num2str(l), 'L_', outLab, '.png']));
    saveas3format(fgrp{l}, GCsummaryDir, ['GCnetGraph_', num2str(l), 'L_', outLab]); %pdf format
end

%% === PLOT STYLE 2: Fixed Arrow Size and Line Width ===
for l = 1:maxLayer
    f = figure;
    g = digraph(plotMat(:,:,l), chanNamesforNet);  % use plotMat here
    p = plot(g, 'layout','layered','Direction','down','Sinks',chanNamesforNet{4},'AssignLayers','alap');
    p.NodeColor = 'k'; p.Marker = 'o'; p.MarkerSize = 8; p.NodeFontSize = 25;
    p.ArrowSize = 25; p.LineWidth = 2;
    p.EdgeAlpha = 0.9;
    edgeColors = zeros(numedges(g),3);
    for e = 1:numedges(g)
        src = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,1}));
        tgt = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,2}));
        if edgeColorMat(src,tgt,l)==1
            edgeColors(e,:) = [1 0 0];
        elseif edgeColorMat(src,tgt,l)==2
            edgeColors(e,:) = [0.5 0.5 0.5];
        end
    end
    p.EdgeColor = edgeColors;
    title(['Granger-Causal Pathway - ', num2str(l), 'L']);
    saveas(f, fullfile(GCsummaryDir, ['GCnetGraph2_', num2str(l), 'L_', outLab, '.png']));
  saveas3format(f, GCsummaryDir, ['GCnetGraph2_', num2str(l), 'L_', outLab]);% pdf format
end

%% === PLOT STYLE 3: Minimal / Binary Edge Rendering ===
for l = 1:maxLayer
    f = figure;
    g = digraph(plotMat(:,:,l), chanNamesforNet);  % use plotMat here
    p = plot(g, 'layout','layered','Direction','down','Sinks',chanNamesforNet{4},'AssignLayers','alap');
    p.NodeColor = 'k'; p.Marker = 'o'; p.MarkerSize = 8; p.NodeFontSize = 25;
    p.ArrowSize = 20; p.LineWidth = 1.5; p.EdgeAlpha = 0.9;
    edgeColors = zeros(numedges(g),3);
    for e = 1:numedges(g)
        src = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,1}));
        tgt = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,2}));
        if edgeColorMat(src,tgt,l)==1
            edgeColors(e,:) = [1 0 0];
        elseif edgeColorMat(src,tgt,l)==2
            edgeColors(e,:) = [0.5 0.5 0.5];
        end
    end
    p.EdgeColor = edgeColors;
    title(['Granger-Causal Pathway - ', num2str(l), 'L']);
    saveas(f, fullfile(GCsummaryDir, ['GCnetGraph3_', num2str(l), 'L_', outLab, '.png']));
   saveas3format(f, GCsummaryDir, ['GCnetGraph3_', num2str(l), 'L_', outLab]);% pdf format
end

%% === PLOT STYLE 4: Red (p<0.05) + Gray (0.05?p<0.10) edges ===
for l = 1:maxLayer
    f = figure;
    g = digraph(plotMat(:,:,l), chanNamesforNet);  % uses both red & gray edge weights
    p = plot(g, 'layout','layered','Direction','down','Sinks',chanNamesforNet{4},'AssignLayers','alap');
    
    % Node appearance
    p.NodeColor = 'k'; p.Marker = 'o'; p.MarkerSize = 8; p.NodeFontSize = 25;
    
    % Edge appearance (arrow size & width based on R²)
    wtVec = g.Edges.Weight;
    p.ArrowSize = 20 + (100-20) * wtVec ./ 100;
    p.LineWidth = 1.5 + (5 - 1.5) * wtVec ./ 100;
    p.EdgeAlpha = 0.9;
    
    % Edge color: red for p<0.05, gray for 0.05?p<0.10
    edgeColors = zeros(numedges(g),3);
    for e = 1:numedges(g)
        src = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,1}));
        tgt = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,2}));
        if edgeColorMat(src,tgt,l) == 1
            edgeColors(e,:) = [1 0 0];         % red
        elseif edgeColorMat(src,tgt,l) == 2
            edgeColors(e,:) = [0.5 0.5 0.5];   % gray
        end
    end
    p.EdgeColor = edgeColors;

    title(['Granger-Causal Pathway - ', num2str(l), 'L']);
    saveas(f, fullfile(GCsummaryDir, ['GCnetGraph4_', num2str(l), 'L_', outLab, '.png']));
    saveas3format(f, GCsummaryDir, ['GCnetGraph4_', num2str(l), 'L_', outLab]);%pdf format
end
%% === PLOT STYLE 5: Label edges with '5%' or '10%' significance thresholds ===
for l = 1:maxLayer
    f = figure;
    g = digraph(plotMat(:,:,l), chanNamesforNet);  % use red + gray edges
    p = plot(g, 'layout','layered','Direction','down','Sinks',chanNamesforNet{4},'AssignLayers','alap');
    % Node appearance
    p.NodeColor = 'k'; p.Marker = 'o'; p.MarkerSize = 8; p.NodeFontSize = 25;
    % Edge appearance
    wtVec = g.Edges.Weight;
    p.ArrowSize = 20 + (100-20) * wtVec ./ 100;
    p.LineWidth = 1.5 + (5 - 1.5) * wtVec ./ 100;
    p.EdgeAlpha = 0.9;
    % Edge color & edge label
    edgeColors = zeros(numedges(g),3);
    edgeLabels = cell(numedges(g), 1);
    for e = 1:numedges(g)
        src = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,1}));
        tgt = find(strcmp(chanNamesforNet, g.Edges.EndNodes{e,2}));
        if edgeColorMat(src,tgt,l) == 1
            edgeColors(e,:) = [1 0 0];       % red
            edgeLabels{e} = '5%';
        elseif edgeColorMat(src,tgt,l) == 2
            edgeColors(e,:) = [0.5 0.5 0.5]; % gray
            edgeLabels{e} = '10%';
        else
            edgeLabels{e} = '';
        end
    end
    p.EdgeColor = edgeColors;
    p.EdgeLabel = edgeLabels;
    p.EdgeFontSize = 15;

    title(['Granger-Causal Pathway - ', num2str(l), 'L']);
    saveas(f, fullfile(GCsummaryDir, ['GCnetGraph5_', num2str(l), 'L_', outLab, '.png']));
    saveas3format(f, GCsummaryDir, ['GCnetGraph5_', num2str(l), 'L_', outLab]);% pdf format
end
%% Save .mat
save(fullfile(GCsummaryDir, ['GCnetworkGraph_Output_', outLab, '.mat']), 'connMat', 'connMatwt', 'maxLayer', 'chanDetailedNames', 'chanNamesforNet');
disp('=== GCsummaryToGCnetGraph4variable (5%/10% with 3 plot styles) DONE ===');

end
