function currSamples = sampleFrameWindowsWithN(movieData, iFrame)
% 
% Modified to add the number of pixels within a window
% by Jungsik Noh
% 6/2016

%% Input

%Make sure the movie has been windowed, and find the desired process.
iWinProc = movieData.getProcessIndex('WindowingProcess');
winProc= movieData.processes_{iWinProc};

iSegProc = winProc.funParams_.SegProcessIndex;   
segProc = movieData.processes_{iSegProc};

%% -------- Init ---------- %%

imSize = movieData.imSize_;
nInput = numel(movieData.channels_);

%% --------- Sampling --------- %%
 
disp(['Starting sampleFrameWindowsWithN...iFrame= ' num2str(iFrame)]);

%iFrame = 1;  
    %Load the windows
    currWin = winProc.loadChannelOutput(iFrame);    
        
    %Load the mask(s) to use first, so we can combine them and use this to 
    %mask every channel.
    currMask = true(imSize);
    %MaskChannelIndex = segProc.funParams_.ChannelIndex;
    %    currMask = currMask & segProc.loadChannelOutput(MaskChannelIndex,iFrame);
        
        %rawIm = movieData.channels_(MaskChannelIndex).loadImage(iFrame);
        %rawImMask = (rawIm > 0);
        %currMask = currMask & rawImMask;       
    
    %Go through each channel and sample it
    stack2sample=zeros([imSize nInput]);
for iChan = 1:nInput
    stack2sample(:,:, iChan) = movieData.channels_(iChan).loadImage(iFrame);
    stack2sample(stack2sample == 0) = NaN;
end
    currSamples = sampleStackWindowsCounts(currWin,stack2sample,currMask);

    
%% Output 

disp(['Finished sampling! iFrame= ' num2str(iFrame)])

end


%% EOF