vars = whos; varsArray = struct2cell(vars); varsArray = varsArray(1,:); varsArray = varsArray(~ismember(varsArray,{'varsArray','origRecon2','origRecon2FVAMin','origRecon2FVAMax'}));
for i=1:length(varsArray)
    clear(varsArray{i});
end
if isunix
    cd /home/ubuntu/MATLAB/PathwayAnalysis
    topDir = '/home/ubuntu/MATLAB/PathwayAnalysis';
else
    cd C:\Users\Yiping' Wang'\Documents\MATLAB\PathwayAnalysis
    topDir = 'C:\Users\Yiping Wang\Documents\MATLAB\PathwayAnalysis';
end
inputDir = [topDir filesep 'input'];
outputDir = [topDir filesep 'output'];

if exist([outputDir filesep 'readPicardData' filesep 'PicardData.mat'],'file')
    load([outputDir filesep 'readPicardData' filesep 'PicardData.mat'],'PicardHeteroplasmies','PicardExpressionData','PicardGeneIDs');
end
if exist([outputDir filesep 'processPicardData' filesep 'PicardDataProcessed.mat'],'file')
    load([outputDir filesep 'processPicardData' filesep 'PicardDataProcessed.mat'],'pValsCorr','pValsT','rhos','labelsCorr','labelsT');
end
uniqueHeteroplasmies = {'0','20','30','50','60','90','100','Rho0'};

if ~exist('origRecon2','var')
    loadModels;
end
if exist([topDir filesep 'FVA.mat'],'file')
    load([topDir filesep 'FVA.mat']);
else
    [origRecon2FVAMin, origRecon2FVAMax, ~, ~] = fluxVariability(origRecon2);
end

configStruct = v2struct();