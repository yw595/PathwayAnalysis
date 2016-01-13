function [sortedNormScores sortedMets sortedMetNames] = RMAPicardDataFunction(geneIDs, pVals, model, sigCutoff, permuteCutoff)
if permuteCutoff==1
    pVals = pVals(randsample(length(pVals),length(pVals)));
end
tempGeneIDs = {}; tempPVals = [];
for j=1:length(geneIDs)
    if any(strcmp(geneIDs{j},model.genes))
	tempGeneIDs{end+1} = geneIDs{j};
        tempPVals(end+1) = pVals(j);
    end
end
geneIDs = tempGeneIDs; pVals = tempPVals;
[sortPVals intersectIdxs] = sort(pVals);
if sigCutoff
    cutoff = FDRCutoff(sortPVals);
else
    cutoff = length(sortPVals);
end
sigGeneIDs = geneIDs(intersectIdxs(1:cutoff)); sigSortPVals = sortPVals(1:cutoff);
if permuteCutoff==2
    sigSortPVals = sigSortPVals(randsample(length(sigSortPVals),length(sigSortPVals)));
end
[sortedNormScores sortedMets sortedMetNames] = runRMA(sigGeneIDs, sigSortPVals, model);

end