function selectExpAllHDData(ensemblIDs, commonNames, configStruct, anOutputDir)

v2struct(configStruct);

count=0;
for i=1:length(uniqueMonths)
    for j=1:length(uniqueQLengths)
        count=count+1;
        xlabels{count} = ['Month' num2str(uniqueMonths(i)) 'QLength' num2str(uniqueQLengths(j))];
    end
end

for l=1:length(uniqueTissues)
    for k=1:length(ensMito)
        disp([uniqueTissues{l} ' ' ensemblIDs{k}]);
        count=0;
        vals = [];
        xvals = [];
        for i=1:length(uniqueMonths)
            for j=1:length(uniqueQLengths)
                count=count+1;
                disp([num2str(i) ' ' num2str(j)])
                [~, intersectIdxs1, intersectIdxs2] = intersect(allGeneIDs, ensemblIDs{k});
                mask = allQLengths == uniqueQLengths(j) & allMonths == uniqueMonths(i) & strcmp(allSeqTypes,'mRNA');
                tempVals = allExpressionData(intersectIdxs1,mask);
                vals = [vals(:); tempVals(:)];
                xvals = [xvals(:); count*ones(length(tempVals(:)),1)];
            end
        end
        makeScatter(xvals, vals, commonNames{k}, xlabels, 'RNASeq Counts', anOutputDir);
    end
end