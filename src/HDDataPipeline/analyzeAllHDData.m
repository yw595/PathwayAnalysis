config;

testT = 0; testCorr = 1;

if(testT)
    pVals = pValsT;
end
if(testCorr)
    pVals = pValsCorr;
end

pVals(isnan(pVals)) = 1;
for i=1:size(pVals,2)
    [sortPValsCol idxs] = sort(pVals(:,i));
    cutoffs(i) = FDRCutoff(sortPValsCol);
    sigEntrezIDs{i} = allGeneIDs(idxs(1:cutoffs(i)));
    sortPValsCol = sortPValsCol(1:cutoffs(i));   
    sortPValsColArray{i} = sortPValsCol;
end

for i=1:size(pVals,2)
    fid = fopen([outputDir filesep 'sigEntrezIDs' num2str(i) '.txt'],'w');
    ithSigEntrezIDs = sigEntrezIDs{i};
    for j=1:length(ithSigEntrezIDs)
        fprintf(fid,'%s\n',ithSigEntrezIDs{j});
    end
    fclose(fid);
end

for k=1:length(uniqueTissues)
    ithSigEntrezIDs = sigEntrezIDs{k};
    if testCorr
        mask = strcmp(allTissues, uniqueTissues{k}) ...
            & allMonths==allMonths(timePoint) & strcmp(allSeqTypes,'mRNA');
    else
        mask = strcmp(allTissues, uniqueTissues{k}) & (allQLengths==20 | allQLengths==175) ...
            & (allMonths==2 | allMonths==10) & strcmp(allSeqTypes,'mRNA');
    end
    [~, intersectIdxs, ~] = intersect(ithLookSigEntrezIDs,allGeneIDs);
    sigEntrezVals{k} = allExpressionData( intersectIdxs,mask );
    sigEntrezQLengths{k} = allQLengths(mask);
end

for i=1:length(lookSigEntrezQLengths)
    lookSigEntrezQLengths{i} = transformQVals(lookSigEntrezQLengths{i});
end

save(outputFile, 'sigEntrezQLengths', 'sigEntrezVals', 'sigEntrezIDs', 'sortPValsColArray');