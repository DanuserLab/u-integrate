function laggedTScolMat = myLagmatrix(TScolMat, lagOrder)
% lagged TS column matrix with intial NaN's imputed by means
%


X = TScolMat;

%%
Xlmat = lagmatrix(X, 1:lagOrder);
Xlmat1 = Xlmat;

% nan impute with mean
for l = 1:lagOrder
    for j = 1:l
        Xlmat1(j, 1+(l-1)*size(X, 2):l*size(X,2)) = mean(X, 1, 'omitnan');
    end
end

laggedTScolMat = Xlmat1;

end
