clear all;
if isunix
    cd /home/ubuntu/MATLAB/PathwayAnalysis
    topDir = '/home/ubuntu/MATLAB/PathwayAnalysis';
else
    cd C:\Users\Yiping' Wang'\Documents\MATLAB\PathwayAnalysis
    topDir = 'C:\Users\Yiping Wang\Documents\MATLAB\PathwayAnalysis';
end
inputDir = [topDir filesep 'input'];
outputDir = [topDir filesep 'output'];

if exist([outputDir filesep 'integrateAllHDData' filesep 'allHDData.mat'],'file')
    load([outputDir filesep 'integrateAllHDData' filesep 'allHDData.mat'],'allObservationIDs','allMouseIDs','allTissues', ...
        'allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');
end
if exist([outputDir filesep 'processAllHDData' filesep 'allHDDataProcessed.mat'],'file')
    load([outputDir filesep 'processAllHDData' filesep 'allHDDataProcessed.mat'],'pValsCorr','pValsT','rhos','labelsCorr','labelsT');
end
uniqueMonths = {'2Months', '6Months', '10Months'};
uniqueQLengths = {'20Q', '80Q', '92Q', '111Q', '140Q', '175Q'};
uniqueTissues = {'cortex','liver','striatum','cerebellum','hippocampus'};
uniqueSexes = {'FEMALE','MALE'};

if ~exist('origRecon2','var')
    loadModels;
end

configStruct = v2struct();