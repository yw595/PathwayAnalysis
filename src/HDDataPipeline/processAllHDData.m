config;
v2struct();
disp('running processAllHDData');
outputDir1 = [outputDir filesep 'processAllHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

pValsCorr = zeros( size(allExpressionData,1),length(uniqueTissues),length(uniqueMonths) );
pValsT = zeros( size(allExpressionData,1),length(uniqueTissues));
rhos = zeros( size(allExpressionData,1),length(uniqueTissues));

labelsT = {uniqueTissues,{'QLength','Month'}};
labelsCorr = {uniqueTissues,uniqueMonths};
for i=1:size(allExpressionData,1)
    i
    if mod(i,1000)==0
        disp(i);
    end
    for k=1:length(uniqueTissues)
        mask = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA');
        dataRow = allExpressionData(i,:);
        mask1 = mask & strcmp(allQLengths,'20Q');
        mask2 = mask & strcmp(allQLengths,'175Q');
        [h pValsT(i,k,1)] = ttest2(dataRow(mask1)', dataRow(mask2)');
        mask1 = mask & strcmp(allMonths,'2Months');
        mask2 = mask & strcmp(allMonths,'10Months');
        [h pValsT(i,k,2)] = ttest2(dataRow(mask1)', dataRow(mask2)');

        for j=1:length(uniqueMonths)
            mask1 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & strcmp(allMonths,uniqueMonths{j});
            [rho pValsCorr(i,k,j)] = corr(cellfun(@(x) str2num(x(1:end-1)), allQLengths(mask1)'), dataRow(mask1)', 'type', 'Spearman');
        end
    end
end

save([outputDir1 filesep 'allHDDataProcessed.mat'],'pValsCorr','pValsT','rhos','labelsT','labelsCorr');