config;
v2struct();
disp('running integrateAllHDData');
outputDir1 = [outputDir filesep 'integrateAllHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

load([outputDir filesep 'readThreeHDData' filesep 'threeHDData.mat']);
load([outputDir filesep 'readCAndHHDData' filesep 'cAndHHDData.mat']);
[cAndHGeneIDs, ~, ~, cAndHExpressionData] = mapGeneIDs([inputDir filesep 'mouseEnsemblToMouseEntrez.csv'], cAndHGeneIDs, 0, 'expressionData', cAndHExpressionData);

[allGeneIDs, intersectIdxs1, intersectIdxs2] = intersect(threeGeneIDs, cAndHGeneIDs);
allExpressionData = [threeExpressionData(intersectIdxs1,:) cAndHExpressionData(intersectIdxs2,:)];
allObservationIDs = [threeObservationIDs' cAndHObservationIDs];
allTissues = [threeTissues' cAndHTissues];
allQLengths = [threeQLengths' cAndHQLengths];
allSexes = [threeSexes' cAndHSexes];
allMonths = [threeMonths' cAndHMonths];
allSeqTypes = [threeSeqTypes' cAndHSeqTypes];

save([outputDir1 filesep 'allHDData.mat'],'allObservationIDs', ...
    'allTissues','allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');