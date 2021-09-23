function MLsummary_FluctuationProfiling(ML, chanName, maxLayer, analNameDesc, outDirName, varargin)
% MLsummary_FluctuationProfiling Collect/Summarize the output from
% phaseDescriptives_OneChan().
%
% Usage:
% MLsummary_FluctuationProfiling(ML, 'Actin', 2, 'phaseDescriptives_0p4mRL5', 'phaseDesc0p4mRL5_Actin')
%
% Updated:
% J Noh, 2018/10/29. Add a 'figFlag' option. 
% J Noh, 2018/02/22. Add a 'timeInterval' option to handle the 'Folding'
% option.
% Jungsik Noh, 2018/01/30


% load ML, example MD
ML.getMovies()
md1 = ML.getMovie(1);
[folderName, cellLab0, ~] = fileparts(md1.outputDirectory_)


%%

ip = inputParser; 
ip.addParameter('lagMax0', nan);
ip.addParameter('timeInterval', md1.timeInterval_);
ip.addParameter('figFlag', 'off');


parse(ip, varargin{:})
p = ip.Results;


%%  Param (need to be converted to arguments)
% basically input arguments

ch0ActmapName = 'Vel';
ch1ActmapName0 = chanName;


outDir = fullfile(ML.outputDirectory_, outDirName);   % analName
if ~isdir(outDir); mkdir(outDir); end

%
ZsamplesName = [ch1ActmapName0, '-ret-ch1Zsamples.mat'];
load(fullfile(md1.outputDirectory_, analNameDesc, ZsamplesName));  % xcorrMat_tmp

size(ch1Zsamples{1}{1})'

if ~isnan(p.lagMax0)
    frSize = 2*p.lagMax0 + 1;
else
    frSize = size(ch1Zsamples{1}{1}, 2);
end

disp(['frSize: ', num2str(frSize)])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

MDs = ML.getMovies();
num = numel(MDs);

%%  Gather .mat files and .png figures
%%  '-ret'

ch1ActmapName = [ch1ActmapName0, '-ret']

MLsummary_phaseDesc_Loop(md1, ch0ActmapName, maxLayer, outDir, ...
                                analNameDesc, MDs, num, ch1ActmapName, frSize, p.timeInterval, p.figFlag)
pause(2)

%%  '-prot'

ch1ActmapName = [ch1ActmapName0, '-prot']

MLsummary_phaseDesc_Loop(md1, ch0ActmapName, maxLayer, outDir, ...
                                analNameDesc, MDs, num, ch1ActmapName, frSize, p.timeInterval, p.figFlag)
pause(2)

%%  '-maxVel'

ch1ActmapName = [ch1ActmapName0, '-maxVel']

MLsummary_phaseDesc_Loop(md1, ch0ActmapName, maxLayer, outDir, ...
                                analNameDesc, MDs, num, ch1ActmapName, frSize, p.timeInterval, p.figFlag)
pause(2)

%%  '-minVel'

ch1ActmapName = [ch1ActmapName0, '-minVel']

MLsummary_phaseDesc_Loop(md1, ch0ActmapName, maxLayer, outDir, ...
                                analNameDesc, MDs, num, ch1ActmapName, frSize, p.timeInterval, p.figFlag)
pause(2)

%%
disp('==== MLsummary_FluctuationProfiling is finished!! ====')
disp('== :) close all ')



end


