function [fileNames] = getFileNames(folder, jobRange)
%getFileNames get a list of job file names in the 'folder' directory
%   Filenames must contain a substring of the format *job[jobNum]*.mat

fileNames = [];
for f_it = 1:numel(jobRange)
    d = dir([folder 'job' num2str(jobNums(f_it)) '*.mat']);


end

