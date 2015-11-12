uniqueTissues = {'cortex','liver','striatum'};
uniqueSexes = {'FEMALE','MALE'};
pValsCorr = zeros( 3,size(allExpressionData,1),length(uniqueTissues) );
pValsT = zeros( size(allExpressionData,1),length(uniqueTissues));
rhos = zeros( size(allExpressionData,1),length(uniqueTissues));

for i=1:size(allExpressionData,1)
    if mod(i,1000)==0
        disp(i);
    end
    count = 0;
    for k=1:length(uniqueTissues)
        count = count+1;
        mask = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA');
        dataRow = allExpressionData(i,:);
        mask1 = mask & allQLengths==20 & allMonths==2;
        mask2 = mask & allQLengths==175 & allMonths==10;
        [h pValsT(i,count)] = ttest2(dataRow(mask1)', dataRow(mask2)');
        
        mask1 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==2;
        [rho pValsCorr(1,i,count)] = corr(allQLengths(mask1), dataRow(mask1)', 'type', 'Spearman');
        mask2 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==6;
        [rho pValsCorr(2,i,count)] = corr(allQLengths(mask2), dataRow(mask2)', 'type', 'Spearman');
        mask3 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==10;
        [rho pValsCorr(3,i,count)] = corr(allQLengths(mask3), dataRow(mask3)', 'type', 'Spearman');
    end
end

save([inputDir filesep 'AllHDDataProcessed.mat'],'uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');