load('HDData.mat');

uniqueTissues = {'cortex','Liver','striatum'};
uniqueSexes = {'F','M'};
pValsCorr = zeros( size(expressionDataHD,1),length(uniqueTissues)*length(uniqueSexes) );
pValsT = zeros( size(expressionDataHD,1),length(uniqueTissues)*length(uniqueSexes) );
rhos = zeros( size(expressionDataHD,1),length(uniqueTissues)*length(uniqueSexes) );
for i=1:size(expressionDataHD,1)
    if mod(i,1000)==0
        disp(i);
    end
    count = 0;
    for k=1:length(uniqueTissues)
        for l=1:length(uniqueSexes)
            count = count+1;
            mask = strcmp(tissues, uniqueTissues{k}) & strcmp(sexes, uniqueSexes{l});
            dataRow = expressionDataHD(i,:); [rho pVal] = ...
                corr(cellfun(@(x) str2num(x(2:end)),QLengths(mask)), dataRow(mask)', 'type', 'Spearman');
            pValsCorr(i,count) = pVal; rhos(i,count) = rho;
            mask1 = mask & strcmp(QLengths, '20');
            mask2 = mask & strcmp(QLengths, '175');
            [h pVal] = ttest2(dataRow(mask1)', dataRow(mask2)'); pValsT(i,count) = pVal;
        end
    end
end

save('HDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');