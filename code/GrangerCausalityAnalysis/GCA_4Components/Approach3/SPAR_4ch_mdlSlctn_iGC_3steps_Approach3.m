function [arMdlAvBICmap, arMdlAvBICcurve, reducedMdlAvBICmap, reducedMdlAvBICcurve, ...
    fullMdlAvBICmap, fullMdlAvBICcurve, Jw0, Jw1, ...
    Pw01, Pw11,  Pw02, Pw12, Pw13,Qw0, modelIndVarsMat, modelCoefMat, winFitsMat, resiMat, yhatFMat, yhatRMat] = ...
    SPAR_4ch_mdlSlctn_iGC_3steps_Approach3(Mapsx, Mapsy, Mapsz, Mapsv, ...
    tlagMax, wlagMax, tlagMaxReg, wlagMaxReg, regPvecMat, ICtype)

wmax = size(Mapsy{2}, 1);
arMdlAvBICmap = nan(wmax, size(regPvecMat, 1)); 
w00 = max(wlagMax, wlagMaxReg);
%NumLags=20;
%% loop over windows for AR model

for w = (1+w00):(wmax-w00)    
    
    fprintf(1, '%g ', w); if (mod(w,50)==0); fprintf('\n'); end
    [yadj, ~, ~,~] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz, Mapsv);
    
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
            
            bhat = X\y; sse = (y-X*bhat)'*(y-X*bhat);  % faster?
            avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
            avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
            avAIC = -2*avlogLhat + 2/N*(size(X,2) + 1);
            
            if strcmp(ICtype, 'AIC')
                tmpBICvec(j) = avAIC;
            elseif strcmp(ICtype, 'BIC')
                tmpBICvec(j) = avBIC;
            else
                error('ICtype must be one of AIC or BIC')
     end
        end
    end
    arMdlAvBICmap(w,:) = tmpBICvec;    
end
%% regression selections 
jmax = tlagMaxReg+1;   
arMdlAvBICcurve = nanmean(arMdlAvBICmap, 1);
[~, i2] = min(arMdlAvBICcurve);
arPvec0 = regPvecMat(i2,:);
Jw0 = arPvec0(1);
Jw1 = arPvec0(2);
%%
fullMdlAvBICmap    = nan(wmax, jmax); 
%%
for w = (1+w00):(wmax-w00)    
    
    fprintf(1, '%g ', w); 
    if (mod(w,50)==0); fprintf('\n'); end
    
    [yadj, xadj, ~, ~] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz, Mapsv);
    
    y = yadj{1};
    N = size(y, 1);
    if all(isnan(y)) || any(isempty(y))
        continue
    end
    %% build AR + Reg. (reduced model)
    yLag = cell(2*1+1+2, 1);
    yLag{1} = myLagmatrix(yadj{1}, Jw0);
    X1 = [yLag{1}];
    
    if wlagMax > 0
        for k = 1:wlagMax
            yLag{2*k}   = myLagmatrix(yadj{2*k},   Jw1);
            yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, Jw1);
            X1 = [X1, yLag{2*k}, yLag{2*k+1}];
        end
    end
    % for layer propagation
    yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, Jw1);
    yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, Jw1);
    X1 = [X1, yLag{end-1}, yLag{end}];
    
    % reduced reg model
    X_ArCtrl = zscore(X1);
    nanIndVars = any(isnan(X1), 1);
    X1sub = X_ArCtrl(:, ~nanIndVars);
   
    %% full regression loop
    tmpBICvec = nan(1, jmax);   % initialize once per w
    
    for j = 1:jmax
        % regPvec
        regPvec = j-1;
        
        % instantaneous GC + lag
        xLag = myLagmatrix(xadj{1}, regPvec);
        X2 = [xadj{1}, xLag];
        
        % handle NaNs
        nanind = any(isnan(X2), 1);
        X_Reg = X2(:, ~nanind);
        
        % combine reduced + regressor
        X = zscore([X1sub, X_Reg]);
        
        if any(any(isnan([y, X])))
            tmpBICvec(j) = nan;
        else
            bhat = X\y; 
            sse = (y - X*bhat)'*(y - X*bhat);
            avlogLhat = -0.5*(1+log(sse/N)) - 0.5*log(2*pi);
            avBIC = -2*avlogLhat + log(N)/N*(size(X,2) + 1);
            avAIC = -2*avlogLhat + 2/N*(size(X,2) + 1);
            
            if strcmp(ICtype, 'AIC')
                tmpBICvec(j) = avAIC;
            elseif strcmp(ICtype, 'BIC')
                tmpBICvec(j) = avBIC;
            else
                error('ICtype must be one of AIC or BIC')
            end
        end
    end % for j
    
    fullMdlAvBICmap(w,:) = tmpBICvec;  
    
end % for w

%% Store optimal lag per window (optional)
 fullMdlAvBICcurve = nanmean(fullMdlAvBICmap, 1);
%[~, optimalIdx] = min(tmpBICvec); % find lag combination with minimum IC
[~, regPvec0] = min(fullMdlAvBICcurve);
regQw0 = regPvec0 - 1;
% find lag combination with minimum IC
 Qw0 = regQw0;
 %% Controls Under the same for loops
%% Initialize map
nJ1 = size(regPvecMat,1); % for Pw01,Pw11
nJ2 = size(regPvecMat,1); % for Pw02,Pw12
nJ3 = tlagMaxReg;         % for Pw13
nComb = nJ1 * nJ2 * nJ3;  % total combinations per window
reducedMdlAvBICmap = nan(wmax, nComb); % store results per window

%% Loop over windows m
for w = (1+w00):(wmax-w00)
    fprintf(1, '%g ', w); if mod(w,50)==0; fprintf('\n'); end
    [yadj, zadj, vadj, xadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz, Mapsv);
    y = yadj{1};
    N = size(y,1);

    if all(isnan(y)) || any(isempty(y))
        continue
    end

    %% Nested loops for control lags
    idx = 0; % linear index for tmpBICvec
    tmpBICvec = nan(1, nComb); % correct size for nested loop
    for j1 = 1:nJ1
        Pw01 = regPvecMat(j1,1);
        Pw11 = regPvecMat(j1,2);

        for j2 = 1:nJ2
            Pw02 = regPvecMat(j2,1);
            Pw12 = regPvecMat(j2,2);

            for j3 = 1:nJ3
                Pw13 = j3;
                idx = idx + 1;

                %% Construct lagged matrices for y
                yLag = cell(2*1+1+2,1);
                yLag{1} = myLagmatrix(yadj{1}, Jw0);
                X1 = [yLag{1}];

                if wlagMax > 0
                    for k = 1:wlagMax
                        yLag{2*k} = myLagmatrix(yadj{2*k}, Jw1);
                        yLag{2*k+1} = myLagmatrix(yadj{2*k+1}, Jw1);
                        X1 = [X1, yLag{2*k}, yLag{2*k+1}];
                    end
                end

                %% Layer propagation
                yLag{1+2*wlagMax+1} = myLagmatrix(yadj{1+2*wlagMax+1}, Jw1);
                yLag{1+2*wlagMax+2} = myLagmatrix(yadj{1+2*wlagMax+2}, Jw1);
                X1 = [X1, yLag{end-1}, yLag{end}];
                %%  Include X
                xLag = cell(2*wlagMaxReg+1+2, 1);
                X1 = [X1, xadj{1}];
                xLag{1} = myLagmatrix(xadj{1}, Qw0);
                X1 = [X1, xLag{1}];

                %% Include controls (Z)
                X1 = [X1, zadj{1}];
                zLag{1} = myLagmatrix(zadj{1}, Pw01);
                X1 = [X1, zLag{1}];

                if wlagMax > 0
                    for k = 1:1
                        X1 = [X1, zadj{2}];
                        zLag{2*k} = myLagmatrix(zadj{2*k}, Pw11);
                        X1 = [X1, zLag{2}];
                        X1 = [X1, zadj{3}];
                        zLag{2*k+1} = myLagmatrix(zadj{2*k+1}, Pw11);
                        X1 = [X1, zLag{2*k+1}];
                    end
                end

                zLag{1+2*wlagMax+1} = myLagmatrix(zadj{1+2*wlagMax+1}, Pw11);
                X1 = [X1, zadj{4}, zLag{4}];
                zLag{1+2*wlagMax+2} = myLagmatrix(zadj{1+2*wlagMax+2}, Pw11);
                X1 = [X1, zadj{5}, zLag{5}];

                %% Include V_t(w,l)
                vLag = cell(2*wlagMax+1+2, 1);
                X1 = [X1, vadj{1}];
                vLag{1} = myLagmatrix(vadj{1}, Pw02);
                X1 = [X1, vLag{1}];

                if wlagMax > 0
                    for k = 1:1
                        X1 = [X1, vadj{2}];
                        vLag{2*k} = myLagmatrix(vadj{2*k}, Pw12);
                        X1 = [X1, vLag{2}];
                        X1 = [X1, vadj{3}];
                        vLag{2*k+1} = myLagmatrix(vadj{2*k+1}, Pw12);
                        X1 = [X1, vLag{2*k+1}];
                    end
                end

                vLag{1+2*wlagMax+1} = myLagmatrix(vadj{1+2*wlagMax+1}, Pw12);
                X1 = [X1, vadj{4}, vLag{4}];
                vLag{1+2*wlagMax+2} = myLagmatrix(vadj{1+2*wlagMax+2}, Pw12);
                X1 = [X1, vadj{5}, vLag{5}];

                %% Controls for mapx
                xcLag = cell(2*wlagMaxReg+1+2, 1);
                if wlagMaxReg > 0
                    for k = 1:1
                        xcLag{2*k} = myLagmatrix(xadj{2*k}, Pw13);
                        X1 = [X1, xadj{2}, xcLag{2}];
                        xcLag{2*k+1} = myLagmatrix(xadj{2*k+1}, Pw13);
                        X1 = [X1, xadj{3}, xcLag{3}];
                    end
                end

                xcLag{1+2*wlagMaxReg+1} = myLagmatrix(xadj{1+2*wlagMaxReg+1}, Pw13);
                X1 = [X1, xadj{4}, xcLag{4}];
                xcLag{1+2*wlagMaxReg+2} = myLagmatrix(xadj{1+2*wlagMaxReg+2}, Pw13);
                X1 = [X1, xadj{5}, xcLag{5}];

                %% Regression + IC
                X_ArCtrl = zscore(X1);
                nanIndVars = any(isnan(X1),1);
                X = X_ArCtrl(:, ~nanIndVars);

                if any(any(isnan([y,X])))
                    tmpBICvec(idx) = nan;
                else
                    bhat = X\y;
                    sse = (y-X*bhat)'*(y-X*bhat);
                    avlogLhat = -1/2*(1+log(sse/N)) - 1/2*log(2*pi);
                    avBIC = -2*avlogLhat + log(N)/N*(size(X,2)+1);
                    avAIC = -2*avlogLhat + 2/N*(size(X,2)+1);

                    if strcmp(ICtype,'AIC')
                        tmpBICvec(idx) = avAIC;
                    elseif strcmp(ICtype,'BIC')
                        tmpBICvec(idx) = avBIC;
                    else
                        error('ICtype must be AIC or BIC');
                    end
                end

            end % j3
        end % j2
    end % j1

    reducedMdlAvBICmap(w,:) = tmpBICvec; % store results for this window
end % w

%% Decode i2 into optimal lags
reducedMdlAvBICcurve = nanmean(reducedMdlAvBICmap, 1);   
[~, i2] = min(reducedMdlAvBICcurve(:));                  
%[j1, j2, j3] = ind2sub([nJ1, nJ2, nJ3], i2); 
[j3, j2, j1] = ind2sub([nJ3, nJ2, nJ1], i2);    
Pw01 = regPvecMat(j1,1);
Pw11 = regPvecMat(j1,2);
Pw02 = regPvecMat(j2,1);
Pw12 = regPvecMat(j2,2);
Pw13 = j3;  % direct index

%% GC F-test & modelIndVarsMat
modelIndVarsMat = nan(wmax, (5)*tlagMax + (5)*(tlagMax+1) + (5)*(tlagMax+1) +(5)*(tlagMaxReg+1));
modelCoefMat = modelIndVarsMat;
winFitsMat = nan(wmax, 15);
resiMat = nan(size(Mapsy{2}));
yhatFMat = nan(size(Mapsy{2}));
yhatRMat = nan(size(Mapsy{2}));

% modelIndVarsFixed
%Pw0 = ctrlPvec0(1); Pw1 = ctrlPvec0(2); 
tmparInd0 = [ones(1, Jw0), zeros(1, tlagMax-Jw0)];
tmparInd1 = [ones(1,Jw1), zeros(1,tlagMax-Jw1)];
tmpIndvec = [tmparInd0, repmat(tmparInd1, 1, 4)];
tmpIndvecReg = [ones(1,Qw0+1), zeros(1,tlagMaxReg-Qw0)];
tmpIndvec2 = [ones(1,Pw01+1), zeros(1,tlagMax-Pw01), ...
    repmat([ones(1,Pw11+1),zeros(1,tlagMax-Pw11)],1,4)]; % variable z

tmpIndvec3 = [ones(1,Pw02+1), zeros(1,tlagMax-Pw02), ...
    repmat([ones(1,Pw12+1),zeros(1,tlagMax-Pw12)],1,4)]; % variable v

tmpIndvec4= repmat([ones(1,Pw13+1),zeros(1,tlagMax-Pw13)],1,4);  % xcLag

modelIndVarsFixed = [tmpIndvec,tmpIndvec2, tmpIndvec3,tmpIndvec4, tmpIndvecReg];

for w = (1+w00):(wmax-w00)
    fprintf(1, '%g ', w); 
    if (mod(w,50)==0); fprintf('\n'); end
    
    [yadj, xadj, zadj, vadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz, Mapsv);
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

    % vlag part of ctr
    vLag = cell(2*wlagMax+1+2, 1);
    X1 = [X1, vadj{1}];
    vLag{1} = myLagmatrix(vadj{1}, tlagMax);
    X1 = [X1, vLag{1}];

    if wlagMax > 0
        for k = 1:1
            X1 = [X1, vadj{2}];
            vLag{2*k} = myLagmatrix(vadj{2*k}, tlagMax);
            X1 = [X1, vLag{2}];
            X1 = [X1, vadj{3}];
            vLag{2*k+1} = myLagmatrix(vadj{2*k+1}, tlagMax);
            X1 = [X1, vLag{2*k+1}];
        end
    end    
    vLag{1+2*wlagMax+1} = myLagmatrix(vadj{1+2*wlagMax+1}, tlagMax);
    X1 = [X1, vadj{4}, vLag{4}];
    vLag{1+2*wlagMax+2} = myLagmatrix(vadj{1+2*wlagMax+2}, tlagMax);
    X1 = [X1, vadj{5}, vLag{5}];

    % control lag part
    xcLag = cell(2*wlagMaxReg+1+2, 1);
    if wlagMaxReg > 0
        for k = 1:1
            xcLag{2*k} = myLagmatrix(xadj{2*k}, tlagMax);
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
    X2 = xadj{1};
    xLag{1} = myLagmatrix(xadj{1}, tlagMaxReg);
    X2 = [X2, xLag{1}]; 

    % full data matrix
    X=[X1, X2];     

    %% modelIndVars2 
nanIndVars = any(isnan(X), 1);
modelIndVars2 = modelIndVarsFixed;
modelIndVars2(nanIndVars) = NaN;
modelIndVarsMat(w, :) = modelIndVars2;
Xchosen = X(:, (modelIndVars2 == 1));

%% nan is not allowed.
if any(any(isnan([y, Xchosen]))) || any(isempty(Xchosen(:)))
    fitOut = table(nan(1, 15));  
    r_f = nan(N,1); yhatF = nan(N,1); yhatR = nan(N,1);
    modelIndVars2 = nan(1, size(X, 2));
    modelIndVarsMat(w, :) = modelIndVars2; 
    continue
end

% Full regression
Xsub = zscore(Xchosen); 
c=corr(Xsub(:, 1:end));
e=eig(c);
es=sort(e, 'descend');
cnum = sqrt(es(1)/es(end));
bhat = Xsub\y;
sse = (y - Xsub*bhat)'*(y - Xsub*bhat);
Rsq0 = 1 - sse/( (y - mean(y))'*(y - mean(y)) );
errVar = sse/(N - size(Xsub, 2));
avlogLhat = -1/2*(1 + log(sse/N)) - 1/2*log(2*pi);
avBIC = -2*avlogLhat + log(N)/N*(size(X, 2) + 1);
avAIC = -2*avlogLhat + 2/N*(size(Xsub, 2) + 1);
r_f = y - Xsub * bhat;
yhatF = Xsub*bhat;

[~, KStestPval_f] = kstest(r_f/std(r_f));
[~, LBpval] = lbqtest(r_f);

% Compute SNR
signalVar = var(yhatF);
noiseVar = var(r_f);
SNR = 10 * log10(signalVar / noiseVar);

% coef vec
modelIndCoef = double(modelIndVars2);
modelIndCoef(modelIndVars2 == 1) = bhat;
% coef out
modelCoefMat(w, :) = modelIndCoef;

%% Reduced model (only AR part + ctrl part)
regPart = modelIndVars2(end - tlagMaxReg - 1 + 1:end);
modelIndVars = modelIndVars2(1:end - tlagMaxReg - 1);
reducedX = X(:, (modelIndVars == 1));

Xsub1 = zscore(reducedX);

bhat2 = Xsub1\y;
sse2 = (y - Xsub1 * bhat2)'*(y - Xsub1 * bhat2);
errVar2 = sse2/(N - size(Xsub1, 2));
yhatR = Xsub1 * bhat2;

Fstat = (sse2 - sse)/sum(regPart, 'omitnan') / errVar;
Fpval = fcdf(Fstat, sum(regPart, 'omitnan'), N - size(Xsub, 2), 'upper');
disp(Fpval);

llratio0 = N/2*log(sse2/sse);
LRTstat = 2 * llratio0;
LRTpval = 1 - gamcdf(LRTstat, sum(regPart, 'omitnan')/2, 2);
partialRsquare = (sse2 - sse)/sse2;
Reg_order= Qw0;
% if X is missing 
if (sum(regPart, 'omitnan') == 0)
    partialRsquare = NaN;LRTstat = NaN;LRTpval = NaN;
end

% output
Reg_order= Qw0;
conditionNum = cnum;
Rsquare_full = Rsq0;
AICoverN = avAIC;
BICoverN = avBIC;
RMSE_full = sqrt(errVar);
resiNormalityKSpval_f = KStestPval_f;
RMSE_reduced = sqrt(errVar2);

fitOut = table(AICoverN, conditionNum, Rsquare_full, RMSE_full, ...
    resiNormalityKSpval_f, RMSE_reduced, Fstat, Fpval, ...
    LRTstat, LRTpval, partialRsquare, LBpval, BICoverN, Reg_order, SNR);

winFitsMat(w, :) = table2array(fitOut);
display(winFitsMat(w, :))
resiMat(w, :) = r_f';
yhatFMat(w, :) = yhatF';
yhatRMat(w, :) = yhatR';
end
end
%% Parse spatially organized TS into reg input

function [yadj, xadj, zadj, vadj] = parseSpatiallyOrganizedTS(w, Mapsy, Mapsx, Mapsz, Mapsv)

wlagMax = 1;
wlagMaxReg = 1;

yadj = cell(1+2*wlagMax+2, 1);  % (w,indL) + left/right/up/down
zadj = cell(1+2*wlagMax+2, 1);
vadj = cell(1+2*wlagMax+2, 1);
xadj = cell(1+2*wlagMaxReg+2, 1);

% extract xadj, yadj

yadj{1} = Mapsy{2}(w, :)';
xadj{1} = Mapsx{2}(w, :)';
zadj{1} = Mapsz{2}(w, :)';% extension uses wlagMax (>= wlagMaxReg)
vadj{1} = Mapsv{2}(w, :)';
%vadj{1} = Mapsv{2}(w, :)';
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

if wlagMax > 0
    for k = 1:wlagMax
        vadj{2*k} = Mapsv{2}(w-k, :)';  % left
        vadj{2*k+1} = Mapsv{2}(w+k, :)'; % right
    end
end
vadj{1+2*wlagMax+1}  = Mapsv{1}(w, :)';  % up
vadj{1+2*wlagMax+2} =  Mapsv{3}(w, :)';  % down

if wlagMaxReg > 0
    for k = 1:wlagMaxReg
        xadj{2*k} = Mapsx{2}(w-k, :)';        %left
        xadj{2*k+1} = Mapsx{2}(w+k, :)';      % right
    end
end
xadj{1+2*wlagMaxReg+1} = Mapsx{1}(w, :)';  % up
xadj{1+2*wlagMaxReg+2} = Mapsx{3}(w, :)';  % down

end


