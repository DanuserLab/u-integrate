function outFilePaths = sampleMovieWindowsWithN(movieData)
%function samplesWithN = sampleMovieWindowsWithN(movieData)
% 
% Modified to add the number of pixels within a window
% by Jungsik Noh
% 6/2016

%% Input

%Make sure the movie has been windowed, and find the desired process.
iWinProc = movieData.getProcessIndex('WindowingProcess');
assert(~isempty(iWinProc),'The movie could not be sampled, because it has not been windowed yet!')
winProc= movieData.processes_{iWinProc};

%Make sure that the windows are okay.
assert(winProc.checkChannelOutput,'The window files for the input movie are not valid!')  

%iSegProc = winProc.funParams_.SegProcessIndex;   
%segProc = movieData.processes_{iSegProc};


%% -------- Init ---------- %%

nFrames = movieData.nFrames_;
%imSize = movieData.imSize_;

nInput = numel(movieData.channels_);

%Initialize sample array
sampledFields = {'avg','std','max','min','med','n', 'effectiveN','sum','nWinMask'};
nFields = numel(sampledFields); 

%% --------- Sampling --------- %%

%tic
samplesCell = cell(nFrames, 1);
parfor t = 1:nFrames
    currSamples = sampleFrameWindowsWithN(movieData, t);
    samplesCell{t} = currSamples;
end
%toc

allSamplesArr = nan([nInput, numel(sampledFields), winProc.nSliceMax_, winProc.nBandMax_, nFrames]);
for i = 1:nInput
for j = 1:nFields    
    for t = 1:nFrames
        currSize = size(samplesCell{t}(i).(sampledFields{1}));
        allSamplesArr(i, j, 1:currSize(1), 1:currSize(2), t) ...
            = samplesCell{t}(i).(sampledFields{j});
    end
end
end

%% Output
    
iWinPackInd = movieData.getPackageIndex('WindowingPackage');

if ~isempty(iWinPackInd)
    tmp = movieData.packages_{iWinPackInd}.outputDirectory_;
else
    tmp = movieData.outputDirectory_;
end
    
outFileDirectory_ = fullfile(tmp, 'window_sampling_WithN');
if ~isdir(outFileDirectory_), mkdir(outFileDirectory_), end

    fname = ['all channels.mat'];
    outFilePaths = fullfile(outFileDirectory_, fname);

allSamplesWithN(nInput, 1) = struct();
for i = 1:nInput
    for j = 1:nFields
            allSamplesWithN(i).(sampledFields{j}) ...
                = reshape(allSamplesArr(i, j, :,:,:), winProc.nSliceMax_, winProc.nBandMax_, nFrames);
    end   
end
save(outFilePaths, 'allSamplesWithN');

    
end


%% EOF