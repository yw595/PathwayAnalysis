load(['output' filesep 'MitoCartaVariablesCorr.mat']);

for i=1:length(origRecon2.genes)
    x = origRecon2.genes{i};
    origRecon2Genes{i} = x(1:regexp(x,'\.')-1);
end
for i=1:length(lookSigEntrezIDsArray)
    i
    lookSigEntrezIDs = lookSigEntrezIDsArray{i};
    sortPValsCol = sortPValsColArray{i};
    rxnList = {};
    pValList = [];
    for j=1:length(origRecon2.rxns)
        [relevantGenes,intersectIdxs,~] = intersect(lookSigEntrezIDs, ...
            origRecon2Genes( origRecon2.rxnGeneMat(j,:)~=0 ));
        if ~isempty(relevantGenes)
            rxnList{end+1} = origRecon2.rxns{j};
            pValList(end+1) = mean(sortPValsCol(intersectIdxs));
        end
    end
    normScoresArray{i} = [];
    if ~isempty(pValList)
        normScoresArray{i} = reporterMets(origRecon2,pValList',1000,1,1,'default',rxnList',0);
    end
    [~, sortIdxs] = sort(normScoresArray{i});
    sortedMetsArray{i} = origRecon2.mets(sortIdxs);
    sortedMetNamesArray{i} = origRecon2.metNames(sortIdxs);
end