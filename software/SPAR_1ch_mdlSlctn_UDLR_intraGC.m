function [arMdlAvBICmap, arMdlAvBICcurve,  ...
    Jw0, Jw1, ...
    modelIndVarsMat, modelCoefMat, winFitsMat] = ...
    SPAR_1ch_mdlSlctn_UDLR_intraGC(Mapsy,  ...
    tlagMax, wlagMax, regPvecMat, ICtype)
% SPAR_1ch_mdlSlctn_UDLR_intraGC() IMPLEMENT the regression fitting and model
% selection to detect intracellular information flows in each window. 
%
% Updates:
% 2019/12/30. Add LRUD gc.
% 2019/12/26. Modify false 0 output of partialRsquare, LRTs when X are
% missing. 
% 2019/06/11.
% J Noh, 2019/04/29. Now upper/lower layers of Chan1 are controlling
% variables. (regPvecMat is not unnecessary).
% J Noh, 2019/04/12. Can select BIC option via 'infoCriterion'.
% J Noh, 2019/04/12. Modified from adjAR_...m
%
% Copyright (C) 2025, Danuser Lab - UTSouthwestern 
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
%w00 = max(wlagMax, wlagMaxReg);
w00 = max(wlagMax);

%% loop over windows for AR model

for w = (1+w00):(wmax-w00)    
    
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    yadj = parseSpatiallyOrganizedTS(w, Mapsy);
    
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

%% 2019/12/30
%% GC F-test & modelIndVarsMat
 
modelIndVarsMat = nan(wmax, (5)*tlagMax );
modelCoefMat = modelIndVarsMat;
% 2019/12/30, LRUD gc
winFitsMat = cell(1,4);
for direction = 1:4
    winFitsMat{direction} = nan(wmax, 13);
end

%resiMat = nan(size(Mapsy{2}));
%yhatFMat = nan(size(Mapsy{2}));
%yhatRMat = nan(size(Mapsy{2}));

%Qw0 = regPvec0;
tmparInd0 = [ones(1, Jw0), zeros(1, tlagMax-Jw0)];
tmparInd1 = [ones(1,Jw1), zeros(1,tlagMax-Jw1)];
tmpIndvec = [tmparInd0, repmat(tmparInd1, 1, 4)];
%tmpIndvec3 = repmat([ones(1,Pw1+1),zeros(1,tlagMax-Pw1)],1,4);  % xcLag
%tmpIndvecReg = [ones(1,Qw0+1), zeros(1,tlagMaxReg-Qw0)];

modelIndVarsFixed = [tmpIndvec];

for w = (1+w00):(wmax-w00)
    %disp(w)
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj] = parseSpatiallyOrganizedTS(w, Mapsy);
    
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

    % 2019/04/29. mapx's l/r/u/d now go into control variables.
    % ...
    % full data matrix
    X = [X1];   % ar part + reg part
    
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
        r_f = nan(N,1); %yhatF = nan(N,1); yhatR = nan(N,1);
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
    %yhatF = Xsub*bhat;
    [~, KStestPval_f] = kstest(r_f/std(r_f));
    [~, LBpval] = lbqtest(r_f);          % default lag = 20
    
    % coef vec
    modelIndCoef = double(modelIndVars2);
    modelIndCoef(modelIndVars2 == 1) = bhat;
    % coef out
    modelCoefMat(w, :) = modelIndCoef;
    
    %% reduced model (SPAR - each of U/D/L/R)
    idvec = zeros(size(modelIndVars2)); 
    
    idL = idvec; idL(tlagMax+1:2*tlagMax) = 1;
    idR = idvec; idR(2*tlagMax+1:3*tlagMax) = 1;
    idU = idvec; idU(3*tlagMax+1:4*tlagMax) = 1;
    idD = idvec; idD(4*tlagMax+1:5*tlagMax) = 1;    
    
    regPartL = modelIndVars2(idL == 1);
    regPartR = modelIndVars2(idR == 1);
    regPartU = modelIndVars2(idU == 1);
    regPartD = modelIndVars2(idD == 1);
    
    %modelIndVarsL = modelIndVars2(idL == 0);
    %modelIndVarsR = modelIndVars2(idR == 0);
    %modelIndVarsU = modelIndVars2(idU == 0);
    %modelIndVarsD = modelIndVars2(idD == 0);
    
    reducedXL = X(:, (modelIndVars2 == 1) & (idL == 0));
    reducedXR = X(:, (modelIndVars2 == 1) & (idR == 0));
    reducedXU = X(:, (modelIndVars2 == 1) & (idU == 0));
    reducedXD = X(:, (modelIndVars2 == 1) & (idD == 0));
    
    Xsub1L = zscore(reducedXL);
    Xsub1R = zscore(reducedXR);
    Xsub1U = zscore(reducedXU);
    Xsub1D = zscore(reducedXD);
    
    fitOutL = reducedModelFitting(y, Xsub, sse, errVar, regPartL, Xsub1L, cnum, Rsq0, ...
        avAIC, avBIC, KStestPval_f, LBpval);
    winFitsMat{1}(w,:) = table2array(fitOutL);
    fitOutR = reducedModelFitting(y, Xsub, sse, errVar, regPartR, Xsub1R, cnum, Rsq0, ...
        avAIC, avBIC, KStestPval_f, LBpval);
    winFitsMat{2}(w,:) = table2array(fitOutR);
    fitOutU = reducedModelFitting(y, Xsub, sse, errVar, regPartU, Xsub1U, cnum, Rsq0, ...
        avAIC, avBIC, KStestPval_f, LBpval);
    winFitsMat{3}(w,:) = table2array(fitOutU);
    fitOutD = reducedModelFitting(y, Xsub, sse, errVar, regPartD, Xsub1D, cnum, Rsq0, ...
        avAIC, avBIC, KStestPval_f, LBpval);
    winFitsMat{4}(w,:) = table2array(fitOutD);    

end


end


%% Parse spatially organized TS into reg input

function [yadj] = parseSpatiallyOrganizedTS(w, Mapsy)
% ParseSpatiallyOrganizedTS
% Mapsy are 3*1 cell array with 3 layer activity maps.

wlagMax = 1;
yadj = cell(1+2*wlagMax+2, 1);  % (w,indL) + left/right/up/down

% extract xadj, yadj
yadj{1} = Mapsy{2}(w , :)';     %

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

end

%% reducedModelFitting

function fitOut = reducedModelFitting(y, Xsub, sse, errVar, regPart, Xsub1, cnum, Rsq0, ...
    avAIC, avBIC, KStestPval_f, LBpval)
% reducedModelFitting

    N = size(y, 1);
    
    bhat2 = Xsub1\y;
    sse2 = (y - Xsub1 * bhat2)'*(y - Xsub1 * bhat2);
    errVar2 = sse2/(N-size(Xsub1,2));
    %yhatR = Xsub1*bhat2;

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
    %
    RMSE_reduced = sqrt(errVar2);

    fitOut = table(AICoverN, conditionNum, Rsquare_full, RMSE_full, ...
        resiNormalityKSpval_f, RMSE_reduced, Fstat, Fpval, ...
        LRTstat, LRTpval, partialRsquare, LBpval, BICoverN);   
   
end

    