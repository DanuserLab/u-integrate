function saveas3format(fig, outputDir, fname)
% saveas3format SAVE a figure in .fig, .eps, and .pdf formats.
%
% Jungsik Noh, 2018/01/23

fullname1 = fullfile(outputDir, [fname, '.fig']);
fullname2 = fullfile(outputDir, [fname, '.eps']);
fullname3 = fullfile(outputDir, [fname, '.pdf']);

saveas(fig, fullname1, 'fig');
saveas(fig, fullname2, 'epsc');
saveas(fig, fullname3, 'pdf');


end
