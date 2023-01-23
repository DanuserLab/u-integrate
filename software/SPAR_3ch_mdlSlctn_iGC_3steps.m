function [arMdlAvBICmap, arMdlAvBICcurve, reducedMdlAvBICmap, reducedMdlAvBICcurve, ...
    fullMdlAvBICmap, fullMdlAvBICcurve, Jw0, Jw1, ...
    Pw0, Pw1, Qw0, modelIndVarsMat, modelCoefMat, winFitsMat, resiMat, yhatFMat, yhatRMat] = ...
    SPAR_3ch_mdlSlctn_iGC_3steps(Mapsx, Mapsy, Mapsz, ...
    tlagMax, wlagMax, tlagMaxReg, wlagMaxReg, regPvecMat, ICtype)
% SPAR_3ch_mdlSlctn_iGC_3steps() IMPLEMENT the regression fitting and model
% selection for the GC inference in each window. 
%
% Updates:
% 2019/12/26. Modify false 0 output of partialRsquare, LRTs when X are
% missing. 
% 2019/06/11.
% J Noh, 2019/04/29. Now upper/lower layers of Chan1 are controlling
% variables. (regPvecMat is not unnecessary).
% J Noh, 2019/04/12. Can select BIC option via 'infoCriterion'.
% J Noh, 2019/04/12. Modified from adjAR_...m
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
 
 
% wlagMax === 1, wlagMaxReg===1
wmax = size(Mapsy{2}, 1);
% for SPAR mdl selection
%regPvecMat1 = regPvecMat(2:end, :);
arMdlAvBICmap = nan(wmax, size(regPvecMat, 1));
w00 = max(wlagMax, wlagMaxReg);

%% loop over windows for AR model

parfor w = (1+w00):(wmax-w00)    
    
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj, ~, ~] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz);
    
    y = yadj{1};
    N = size(y, 1);
    if all(isnan(y)) || any(isempty(y))
        continue
    end
    %% SPAR (using BIC)
    
    tmpBICvec = nan(1, size(regPvecMat, 1));
    for j = 1:size(regPvecMat, 1)
        %% parse ar+ctrl part
        j0 = regPvecMat(j, 1);
        j1 = regPvecMat(j, 2);
        
        yLag = cell(2*1+1+2, 1);
        yLag{1} = myLagmatrix(yadj{1}, j0);
        X1 = [yLag{1}];
        if wlagMax > 0
            for k = 1:wlagMax
                yLag{2*k} = myLagmatrix(yadj{2*k}, j1);
                yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, j1);
                X1 = [X1, yLag{2*k}, yLag{2*k+1}];
            end
        end
        % for layer propagation
        yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, j1);
        yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, j1);
        X1 = [X1, yLag{end-1}, yLag{end}];
         
        %% SPAR model
        X_Ar = zscore(X1);     % ar + contrl
        % for the case with nan in U/D
        nanIndVars = any(isnan(X1), 1);
        X = X_Ar(:, ~nanIndVars);
        
        % OLS estimate
        % nan is not allowed.
        if any(any(isnan([y, X])))
            tmpBICvec(j) = nan;
            %bhat = nan;
        else
            
            %[~,~,~,~,stats] = regress(y, X);
            %sse = stats(4)*(N-size(X,2));
            %if ~isempty(X); cnum = sqrt(cond(corr(X))); end
            bhat = X\y; sse = (y-X*bhat)'*(y-X*bhat);  % faster?
            avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
            avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
            avAIC = -2*avlogLhat + 2/N*(size(X,2) + 1);
            %avAICc = -2*avlogLhat + 2/N*(size(X,2) + 1) * N/(N- size(X,2)-1 -1);
            
            %fitOut = avAIC;     %, cnum, Rsq0, sqrt(errVar)];  % avBIC, condition num, R2, sigma
            if strcmp(ICtype, 'AIC')
                tmpBICvec(j) = avAIC;
            elseif strcmp(ICtype, 'BIC')
                tmpBICvec(j) = avBIC;
            else
                error('ICtype must be one of AIC or BIC')
            end
        end
        %tmpBICvec(j) = fitOut1;
    end
    arMdlAvBICmap(w,:) = tmpBICvec;    

end

%% SPAR selection

arMdlAvBICcurve = nanmean(arMdlAvBICmap, 1);
%figure, plot(arMdlAvBICmap')
% figure, plot(arMdlAvBICcurve)
[~, i2] = min(arMdlAvBICcurve);
arPvec0 = regPvecMat(i2,:);
% ar chosen order
Jw0 = arPvec0(1);
Jw1 = arPvec0(2);


%% loop over windows for reducedMdlAvBICmap

reducedMdlAvBICmap = nan(wmax, size(regPvecMat, 1));

parfor w = (1+w00):(wmax-w00)    
    %disp(w)
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj, xadj, zadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz);
    
    y = yadj{1};
    N = size(y, 1);
    %%
    if all(isnan(y)) || any(isempty(y))
        %    minAvAIC2 = nan;
        continue
    end
    
    %% adj AR for best reduced models (using BIC)
    
    tmpBICvec = nan(1, size(regPvecMat, 1));
    for j = 1:size(regPvecMat, 1)
        %% parse ar+ctrl part
        Pw0 = regPvecMat(j, 1);
        Pw1 = regPvecMat(j, 2);
        
        yLag = cell(2*1+1+2, 1);
        yLag{1} = myLagmatrix(yadj{1}, Jw0);
        X1 = [yLag{1}];
        if wlagMax > 0
            for k = 1:wlagMax
                yLag{2*k} = myLagmatrix(yadj{2*k}, Jw1);
                yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, Jw1);
                X1 = [X1, yLag{2*k}, yLag{2*k+1}];
            end
        end
        % for layer propagation
        yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, Jw1);
        yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, Jw1);
        X1 = [X1, yLag{end-1}, yLag{end}];
        
        zLag = cell(2*wlagMax+1+2, 1);
        % include Z_t(w,l) and all adjacent signals at t
        X1 = [X1, zadj{1}];
        zLag{1} = myLagmatrix(zadj{1}, Pw0);
        X1 = [X1, zLag{1}];
        if wlagMax > 0
            for k = 1:1
                X1 = [X1, zadj{2}];
                zLag{2*k} = myLagmatrix(zadj{2*k}, Pw1);
                X1 = [X1, zLag{2}];
                X1 = [X1, zadj{3}];
                zLag{2*k+1} = myLagmatrix(zadj{2*k+1}, Pw1);
                X1 = [X1, zLag{2*k+1}];
            end
        end
        zLag{1+2*wlagMax+1} = myLagmatrix(zadj{1+2*wlagMax+1}, Pw1);
        X1 = [X1, zadj{4}, zLag{4}];
        zLag{1+2*wlagMax+2} = myLagmatrix(zadj{1+2*wlagMax+2}, Pw1);
        X1 = [X1, zadj{5}, zLag{5}];
        
        % 2019/04/29. mapx's l/r/u/d now go into control variables.
        xcLag = cell(2*wlagMaxReg+1+2, 1);
        %xcLag{1} = myLagmatrix(zadj{1}, tlagMax);
        %X1 = [X1, xcLag{1}];
        if wlagMaxReg > 0
            for k = 1:1
                xcLag{2*k} = myLagmatrix(xadj{2*k}, Pw1);     % tlagMax because it is controling
                X1 = [X1, xadj{2}, xcLag{2}];
                xcLag{2*k+1} = myLagmatrix(xadj{2*k+1}, Pw1);
                X1 = [X1, xadj{3}, xcLag{3}];
            end
        end
        xcLag{1+2*wlagMaxReg+1} = myLagmatrix(xadj{1+2*wlagMaxReg+1}, Pw1);
        X1 = [X1, xadj{4}, xcLag{4}];
        xcLag{1+2*wlagMaxReg+2} = myLagmatrix(xadj{1+2*wlagMaxReg+2}, Pw1);
        X1 = [X1, xadj{5}, xcLag{5}];
        
        %% reduced reg model
        X_ArCtrl = zscore(X1);     % ar + contrl
        % for the case with nan in U/D
        nanIndVars = any(isnan(X1), 1);
        X = X_ArCtrl(:, ~nanIndVars);
        
        % OLS estimate
        % nan is not allowed.
        if any(any(isnan([y, X])))
            tmpBICvec(j) = nan;
            %bhat = nan;
        else
            
            %[~,~,~,~,stats] = regress(y, X);
            %sse = stats(4)*(N-size(X,2));
            %if ~isempty(X); cnum = sqrt(cond(corr(X))); end
            bhat = X\y; sse = (y-X*bhat)'*(y-X*bhat);  % faster?
            avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
            avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
            avAIC = -2*avlogLhat + 2/N*(size(X,2) + 1);
            %avAICc = -2*avlogLhat + 2/N*(size(X,2) + 1) * N/(N- size(X,2)-1 -1);
            
            %fitOut = avAIC;     %, cnum, Rsq0, sqrt(errVar)];  % avBIC, condition num, R2, sigma
            if strcmp(ICtype, 'AIC')
                tmpBICvec(j) = avAIC;
            elseif strcmp(ICtype, 'BIC')
                tmpBICvec(j) = avBIC;
            else
                error('ICtype must be one of AIC or BIC')
            end
        end
        %tmpBICvec(j) = fitOut1;
    end
    reducedMdlAvBICmap(w,:) = tmpBICvec;
end

%% reduced mdl selection

reducedMdlAvBICcurve = nanmean(reducedMdlAvBICmap, 1);
%figure, imagesc(reducedMdlAvBICmap), colormap(jet)
%figure, plot(reducedMdlAvBICcurve);
%figure, plot(reducedMdlAvBICmap')
[~, i2] = min(reducedMdlAvBICcurve);
ctrlPvec0 = regPvecMat(i2,:);

%% loop over windows for full model

%regAvAICvec = nan(tlagMaxReg+1, 1);
fullMdlAvBICmap = nan(wmax, tlagMaxReg+1);
% parse ar+ctrl part
Pw0 = ctrlPvec0(1);
Pw1 = ctrlPvec0(2);

jmax = tlagMaxReg+1;

parfor w = (1+w00):(wmax-w00)
    %disp(w)
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj, xadj, zadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz);
    
    y = yadj{1};
    N = size(y, 1);
    
    % parse chosen ar+ctrl model with Pw0, Pw1
    yLag = cell(2*1+1+2, 1);
    yLag{1} = myLagmatrix(yadj{1}, Jw0);
    X1 = [yLag{1}];
    if wlagMax > 0
        for k = 1:wlagMax
            yLag{2*k} = myLagmatrix(yadj{2*k}, Jw1);
            yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, Jw1);
            X1 = [X1, yLag{2*k}, yLag{2*k+1}];
        end
    end
    % for layer propagation
    yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, Jw1);
    yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, Jw1);
    X1 = [X1, yLag{end-1}, yLag{end}];
    
    zLag = cell(2*wlagMax+1+2, 1);
    % include Z_t(w,l) and all adjacent signals at t
    X1 = [X1, zadj{1}];
    zLag{1} = myLagmatrix(zadj{1}, Pw0);
    X1 = [X1, zLag{1}];
    if wlagMax > 0
        for k = 1:1
            X1 = [X1, zadj{2}];
            zLag{2*k} = myLagmatrix(zadj{2*k}, Pw1);
            X1 = [X1, zLag{2}];
            X1 = [X1, zadj{3}];
            zLag{2*k+1} = myLagmatrix(zadj{2*k+1}, Pw1);
            X1 = [X1, zLag{2*k+1}];
        end
    end
    zLag{1+2*wlagMax+1} = myLagmatrix(zadj{1+2*wlagMax+1}, Pw1);
    X1 = [X1, zadj{4}, zLag{4}];
    zLag{1+2*wlagMax+2} = myLagmatrix(zadj{1+2*wlagMax+2}, Pw1);
    X1 = [X1, zadj{5}, zLag{5}];
    
    % 2019/04/29. mapx's l/r/u/d now go into control variables.
    xcLag = cell(2*wlagMaxReg+1+2, 1);
    %xcLag{1} = myLagmatrix(zadj{1}, tlagMax);
    %X1 = [X1, xcLag{1}];
    if wlagMaxReg > 0
        for k = 1:1
            xcLag{2*k} = myLagmatrix(xadj{2*k}, Pw1);     % tlagMax because it is controling
            X1 = [X1, xadj{2}, xcLag{2}];
            xcLag{2*k+1} = myLagmatrix(xadj{2*k+1}, Pw1);
            X1 = [X1, xadj{3}, xcLag{3}];
        end
    end
    xcLag{1+2*wlagMaxReg+1} = myLagmatrix(xadj{1+2*wlagMaxReg+1}, Pw1);
    X1 = [X1, xadj{4}, xcLag{4}];
    xcLag{1+2*wlagMaxReg+2} = myLagmatrix(xadj{1+2*wlagMaxReg+2}, Pw1);
    X1 = [X1, xadj{5}, xcLag{5}];
    
    %% reduced reg model
    X_ArCtrl = zscore(X1);     % ar + contrl
    % for the case with nan in U/D
    nanIndVars = any(isnan(X1), 1);
    X1sub = X_ArCtrl(:, ~nanIndVars);
    
    %% full reg model 
    for j = 1:jmax
        % regPvec
        regPvec = j-1;
        xLag = cell(1, 1);
        % instantaneous GC, log0 included
        X2 = xadj{1};
        xLag{1} = myLagmatrix(xadj{1}, regPvec);
        X2 = [X2, xLag{1}];
        
        % for the case with nan in U/D (unnecessary to remove)
        nanind = any(isnan(X2), 1);
        X_Reg = X2(:, ~nanind);
        
        X = zscore([X1sub, X_Reg]);   % ar part (optmal) + reg part
        % OLS estimate
        % nan is not allowed.
        if any(any(isnan([y, X])))
            fullMdlAvBICmap(w, j) = nan;
            %bhat = nan;
        else
            %[~,~,~,~,stats] = regress(y, X);
            %sse = stats(4)*(N-size(X,2));
            %if ~isempty(X); cnum = sqrt(cond(corr(X))); end
            bhat = X\y; sse = (y-X*bhat)'*(y-X*bhat);  % faster?
            avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
            avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
            avAIC = -2*avlogLhat + 2/N*(size(X,2) + 1);
            %avAICc = -2*avlogLhat + 2/N*(size(X,2) + 1) * N/(N- size(X,2)-1 -1);
            
            %fitOut = avAIC;     %, cnum, Rsq0, sqrt(errVar)];  % avBIC, condition num, R2, sigma
            if strcmp(ICtype, 'AIC')
                fullMdlAvBICmap(w, j) = avAIC;
            elseif strcmp(ICtype, 'BIC')
                fullMdlAvBICmap(w, j) = avBIC;
            else
                error('ICtype must be one of AIC or BIC')
            end
        end
        %fullMdlAvBICmap(w, j) = fitOut;
    end    
end

%% order selection

fullMdlAvBICcurve = nanmean(fullMdlAvBICmap, 1); 
[minAvIC_full, i2] = min(fullMdlAvBICcurve);
regPvec0 = i2 - 1;

%figure, imagesc(fullMdlAvBICmap), colormap(jet)
%figure, plot(nanmean(fullMdlAvBICmap, 1))
%figure, plot(fullMdlAvBICmap')


%% GC F-test & modelIndVarsMat
 
modelIndVarsMat = nan(wmax, (5)*tlagMax + (5)*(tlagMax+1) + (5)*(tlagMaxReg+1));
modelCoefMat = modelIndVarsMat;
winFitsMat = nan(wmax, 13);
resiMat = nan(size(Mapsy{2}));
yhatFMat = nan(size(Mapsy{2}));
yhatRMat = nan(size(Mapsy{2}));

% modelIndVarsFixed
%Pw0 = ctrlPvec0(1); Pw1 = ctrlPvec0(2);
Qw0 = regPvec0;
tmparInd0 = [ones(1, Jw0), zeros(1, tlagMax-Jw0)];
tmparInd1 = [ones(1,Jw1), zeros(1,tlagMax-Jw1)];
tmpIndvec = [tmparInd0, repmat(tmparInd1, 1, 4)];
tmpIndvec2 = [ones(1,Pw0+1), zeros(1,tlagMax-Pw0), ...
    repmat([ones(1,Pw1+1),zeros(1,tlagMax-Pw1)],1,4)];
tmpIndvec3 = repmat([ones(1,Pw1+1),zeros(1,tlagMax-Pw1)],1,4);  % xcLag
tmpIndvecReg = [ones(1,Qw0+1), zeros(1,tlagMaxReg-Qw0)];

modelIndVarsFixed = [tmpIndvec, tmpIndvec2, tmpIndvec3, tmpIndvecReg];

for w = (1+w00):(wmax-w00)
    %disp(w)
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj, xadj, zadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz);
    
    y = yadj{1};
    N = size(y, 1);
    
    %% full model data matrix X
    yLag = cell(2*wlagMax+1+2, 1);
    yLag{1} = myLagmatrix(yadj{1}, tlagMax);
    X1 = [yLag{1}];
    if wlagMax > 0
        for k = 1:wlagMax
            yLag{2*k} = myLagmatrix(yadj{2*k}, tlagMax);
            yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, tlagMax);
            X1 = [X1, yLag{2*k}, yLag{2*k+1}];
        end
    end
    % for layer propagation
    yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, tlagMax);
    yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, tlagMax);
    X1 = [X1, yLag{end-1}, yLag{end}];

    % ctrling part
    zLag = cell(2*wlagMax+1+2, 1);
    % include Z_t(w,l) and all adjacent signals at t
    X1 = [X1, zadj{1}];
    zLag{1} = myLagmatrix(zadj{1}, tlagMax);
    X1 = [X1, zLag{1}];
    if wlagMax > 0
        for k = 1:1
            X1 = [X1, zadj{2}];
            zLag{2*k} = myLagmatrix(zadj{2*k}, tlagMax);
            X1 = [X1, zLag{2}];
            X1 = [X1, zadj{3}];
            zLag{2*k+1} = myLagmatrix(zadj{2*k+1}, tlagMax);
            X1 = [X1, zLag{2*k+1}];
        end
    end
    zLag{1+2*wlagMax+1} = myLagmatrix(zadj{1+2*wlagMax+1}, tlagMax);
    X1 = [X1, zadj{4}, zLag{4}];
    zLag{1+2*wlagMax+2} = myLagmatrix(zadj{1+2*wlagMax+2}, tlagMax);
    X1 = [X1, zadj{5}, zLag{5}];

    % 2019/04/29. mapx's l/r/u/d now go into control variables.
    xcLag = cell(2*wlagMaxReg+1+2, 1);
    %xcLag{1} = myLagmatrix(zadj{1}, tlagMax);
    %X1 = [X1, xcLag{1}];
    if wlagMaxReg > 0
        for k = 1:1
            xcLag{2*k} = myLagmatrix(xadj{2*k}, tlagMax);     % tlagMax because it is controling
            X1 = [X1, xadj{2}, xcLag{2}];
            xcLag{2*k+1} = myLagmatrix(xadj{2*k+1}, tlagMax);
            X1 = [X1, xadj{3}, xcLag{2*k+1}];
        end
    end
    xcLag{1+2*wlagMaxReg+1} = myLagmatrix(xadj{1+2*wlagMaxReg+1}, tlagMax);
    X1 = [X1, xadj{4}, xcLag{4}];
    xcLag{1+2*wlagMaxReg+2} = myLagmatrix(xadj{1+2*wlagMaxReg+2}, tlagMax);
    X1 = [X1, xadj{5}, xcLag{5}];

    % reg part
    xLag = cell(1, 1);
    % include X_t(w,l)
        X2 = xadj{1};
    xLag{1} = myLagmatrix(xadj{1}, tlagMaxReg);
    X2 = [X2, xLag{1}];    
    % full data matrix
    X = [X1, X2];   % ar part + reg part
    
    %% modelIndVars2    
    %modelIndVars2 = nan(1, size(X, 2));     % tlagMax * 14

    % for the case with nan in U/D
    nanIndVars = any(isnan(X), 1);

    %  
    %modelIndVars2 = ~nanIndVars & modelIndVarsFixed;    
    modelIndVars2 = modelIndVarsFixed;
    modelIndVars2(nanIndVars) = NaN;
    modelIndVarsMat(w, :) = modelIndVars2;
    
    Xchosen = X(:, (modelIndVars2 == 1));
    
    %% nan is not allowed.
    if any(any(isnan([y, Xchosen]))) || any(isempty(Xchosen(:)))
        fitOut = table(nan(1, 13));
        r_f = nan(N,1); yhatF = nan(N,1); yhatR = nan(N,1);
        modelIndVars2 = nan(1, size(X, 2));  
        modelIndVarsMat(w, :) = modelIndVars2;        
        continue
    end
    
    %% full reg
    Xsub = zscore(Xchosen); 
 
    c=corr(Xsub(:, 1:end));
    e=eig(c);     %sum(e)
    es=sort(e, 'descend');
    cnum = sqrt(es(1)/es(end));

    %[~,~,~,~,stats] = regress(y, X);
    %sse = stats(4)*(N-size(X,2));
    bhat = Xsub\y; sse = (y-Xsub*bhat)'*(y-Xsub*bhat);  % faster?
    Rsq0 = 1 - sse/( (y-mean(y))'*(y-mean(y)) );
    errVar = sse/(N-size(Xsub,2));
    avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
    avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
    avAIC = -2*avlogLhat + 2/N*(size(Xsub,2) + 1);
    %avAICc = -2*avlogLhat + 2/N*(size(X,2) + 1) * N/(N- size(X,2)-1 -1);

    r_f = y - Xsub * bhat;
    yhatF = Xsub*bhat;
    [~, KStestPval_f] = kstest(r_f/std(r_f));
    [~, LBpval] = lbqtest(r_f);          % default lag = 20
    
    % coef vec
    modelIndCoef = double(modelIndVars2);
    modelIndCoef(modelIndVars2 == 1) = bhat;
    % coef out
    modelCoefMat(w, :) = modelIndCoef;
    
    %% reduced model (only AR part + ctrl part)
    regPart = modelIndVars2(end-tlagMaxReg-1+1:end);
    modelIndVars = modelIndVars2(1:end-tlagMaxReg-1);
    reducedX = X(:, (modelIndVars == 1));
    Xsub1 = zscore(reducedX);
     
    bhat2 = Xsub1\y;
    sse2 = (y - Xsub1 * bhat2)'*(y - Xsub1 * bhat2);
    errVar2 = sse2/(N-size(Xsub1,2));
    yhatR = Xsub1*bhat2;

    Fstat = (sse2 - sse)/sum(regPart, 'omitnan') / errVar;
    Fpval = fcdf(Fstat, sum(regPart, 'omitnan'), N-size(Xsub,2), 'upper');

    llratio0 = N/2*log(sse2/sse);                  % N*avlogLhat - N*avlogLhat_re
    LRTstat = 2* llratio0;
    LRTpval = 1 - gamcdf(LRTstat, sum(regPart, 'omitnan')/2, 2);    % chi2cdf has a bug!!!! gamcdf is the standard format.
    partialRsquare = (sse2 - sse)/sse2;
    % 2019/12/26. when X are missing. 
    if (sum(regPart, 'omitnan') == 0); 
        partialRsquare = NaN; LRTstat = NaN; LRTpval = NaN;
    end

    % output
    conditionNum = cnum;
    Rsquare_full = Rsq0;
    AICoverN = avAIC;
    BICoverN = avBIC;
    RMSE_full = sqrt(errVar);
    resiNormalityKSpval_f = KStestPval_f;
    RMSE_reduced = sqrt(errVar2);

    fitOut = table(AICoverN, conditionNum, Rsquare_full, RMSE_full, ...
        resiNormalityKSpval_f, RMSE_reduced, Fstat, Fpval, ...
        LRTstat, LRTpval, partialRsquare, LBpval, BICoverN);   
    
    winFitsMat(w,:) = table2array(fitOut);
    resiMat(w,:) = r_f';
    yhatFMat(w,:) = yhatF';
    yhatRMat(w,:) = yhatR';
end

end


%% Parse spatially organized TS into reg input

function [yadj, xadj, zadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz)
% ParseSpatiallyOrganizedTS
% Mapsy are 3*1 cell array with 3 layer activity maps.

wlagMax = 1;
wlagMaxReg = 1;

yadj = cell(1+2*wlagMax+2, 1);  % (w,indL) + left/right/up/down
zadj = cell(1+2*wlagMax+2, 1);
xadj = cell(1+2*wlagMaxReg+2, 1);

% extract xadj, yadj

yadj{1} = Mapsy{2}(w , :)';     %
xadj{1} = Mapsx{2}(w , :)';     % extension uses wlagMax (>= wlagMaxReg)
zadj{1} = Mapsz{2}(w , :)';

% l=1                        y{6}
% l=2              y{4} y{2} y{1} y{3} y{5}
% l=3                        y{7}
if wlagMax > 0
    for k = 1:wlagMax
        yadj{2*k} = Mapsy{2}(w-k, :)';    % left
        yadj{2*k+1} = Mapsy{2}(w+k, :)';  % right
    end
end
yadj{1+2*wlagMax+1} = Mapsy{1}(w, :)';  % up
yadj{1+2*wlagMax+2} = Mapsy{3}(w, :)';  % down


if wlagMax > 0
    for k = 1:wlagMax
        zadj{2*k} = Mapsz{2}(w-k, :)';    % left
        zadj{2*k+1} = Mapsz{2}(w+k, :)';  % right
    end
end
zadj{1+2*wlagMax+1} = Mapsz{1}(w, :)';  % up
zadj{1+2*wlagMax+2} = Mapsz{3}(w, :)';  % down

if wlagMaxReg > 0
    for k = 1:wlagMaxReg
        xadj{2*k} = Mapsx{2}(w-k, :)';        %
        xadj{2*k+1} = Mapsx{2}(w+k, :)';
    end
end
xadj{1+2*wlagMaxReg+1} = Mapsx{1}(w, :)';  % up
xadj{1+2*wlagMaxReg+2} = Mapsx{3}(w, :)';  % down

end
