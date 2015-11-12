clear all;
cd C:\Users\Yiping' Wang'\Documents\MATLAB\PathwayAnalysis
topDir = 'C:\Users\Yiping Wang\Documents\MATLAB\PathwayAnalysis';
inputDir = [topDir filesep 'input'];

load('AllHDData.mat','allObservationIDs','allMouseIDs','allTissues', ...
    'allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');
load('AllHDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');
uniqueMonths = [2, 6, 10];
uniqueQLengths = [20, 80, 92, 111, 140, 175];

configStruct = v2struct();