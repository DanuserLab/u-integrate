function MLsummary_SNRoverSmoothingParams(ML, iChan, chanName, SNRdirName, varargin)
% MLsummary_SNRoverSmoothingParams() SUMMARIZE SNR curves for MDs obtained
% by MD_SNRoverSmoothingParams_EWMA().
%
%
% Updates:
% J Noh, 2019/04/18. Remove Gaussian filter 3by3 smoothing case. 
% J Noh, 2019/04/16.

%% parameter 

ip = inputParser;
%ip.addParameter('ratioThreshold', 0.5);
ip.addParameter('outDirName', SNRdirName);
ip.addParameter('figFlag', 'on');

parse(ip, varargin{:});
p = ip.Results;

% outDir
outDir = fullfile(ML.outputDirectory_, p.outDirName);
if ~isdir(outDir); mkdir(outDir); end

%% load ML, example MD
%ML.getMovies()
%md1 = ML.getMovie(1);
S = load(ML.movieDataFile_{1});
md1 = S.MD;
[~, cellLab0, ~] = fileparts(md1.outputDirectory_);

disp(['The label of the 1st movie will be the folder name for movieData.mat: ', cellLab0])

fname0 = ['Chan', num2str(iChan)];
disp(fname0)

%%
%MDs = ML.getMovies();
%num = numel(MDs);

num = numel(ML.movieDataFile_);
MDs = cell(1, num);
for i = 1:num
    S = load(ML.movieDataFile_{i});
    MDs{i} = S.MD;
end

cellLabels = cell(num, 1);

fname = fullfile(MDs{1}.outputDirectory_, SNRdirName, [fname0, '_medSNR.mat']);
S = load(fname);
medsnr = S.medSNR;
maxLayer = numel(medsnr);
numSms = numel(medsnr{1});
% penalization on roughness
xtickSmPars = [1 - S.smParVec];
smParVec = S.smParVec;

medSNRArr = nan(numSms, num, maxLayer);

for i=1:num
    %
    [~, cellLab0, ~] = fileparts(MDs{i}.outputDirectory_);
    cellLabels{i} = [cellLab0(1:end)]; 
    %
    fname = fullfile(MDs{i}.outputDirectory_, SNRdirName, [fname0, '_medSNR.mat']);
    S = load(fname);
    
    for indL = 1:maxLayer
        medSNRArr(:, i, indL) = S.medSNR{indL};
    end
end

%
save(fullfile(outDir, [fname0, '_medSNRArr.mat']), 'medSNRArr') 

%%  plot

f1 = cell(maxLayer,1);
xaxis0 = reshape(1:numel(xtickSmPars), [], 1);
ymax = max(medSNRArr(:));
ymin = min(medSNRArr(:));

for indL = 1:maxLayer
    
    f1{indL} = figure('Visible', p.figFlag);  
    perCellLines = squeeze(medSNRArr(:, :, indL));
    pp = plot( log10(perCellLines), 'LineWidth', 2);
    ylim([log10(ymin), log10(ymax)])
    
    refline([0,0])
    % snr = 2 or 1/2
    h1 = refline([0, log10(2)]); h1.LineStyle = '--';
    h2 = refline([0, log10(1/2)]); h2.LineStyle = '--';     
    
    set(gca, 'Xtick', 1:numel(xtickSmPars))
    set(gca, 'XTickLabel', xtickSmPars)
    set(gca, 'XTickLabelRotation', 45)
    set(gca, 'FontSize', 12)    
    grid on
    
    xlabel('Smoothness')
    ylabel('log10(SNR) (std(Sig)/std(err))')
    title(['medSNR, ', chanName, '-', num2str(indL), 'L',])

    lgd = legend(cellLabels{1:num}, 'Location', 'eastoutside');
    %lgd.FontSize = 6;
end

%%
for indL = 1:maxLayer
    saveas(f1{indL}, fullfile(outDir, [fname0, '_medSNR_', num2str(indL), 'L.png']), 'png')
end

%%



end
