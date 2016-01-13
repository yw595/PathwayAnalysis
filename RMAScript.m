load(['output' filesep 'MATFilesMitoCartaVariablesCorr.mat']);

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
%     normScoresArray{i} = [];
%     if ~isempty(pValList)
%         normScoresArray{i} = reporterMets(origRecon2,pValList',1000,1,1,'default',rxnList',0);
%     end
%     [~, sortIdxs] = sort(normScoresArray{i});
%     sortedMetsArray{i} = origRecon2.mets(sortIdxs);
%     sortedMetNamesArray{i} = origRecon2.metNames(sortIdxs);
    
    fid = fopen(['PVals' num2str(i) '.txt'],'w');
    for j=1:length(lookSigEntrezIDs)
        fprintf(fid,'%s\t%f\n',lookSigEntrezIDs{j},sortPValsCol(j));
    end
    fclose(fid);
    
    fid = fopen(['Interaction' num2str(i) '.txt'],'w');
    for j=1:length(origRecon2.mets)
        interactingRxns = origRecon2.rxns(origRecon2.S(j,:)~=0);
        for k=1:length(interactingRxns)
            fprintf(fid,'%s\t%s\t%s\n',origRecon2.mets{j},'mr',interactingRxns{k});
        end
    end
    fclose(fid);
    
    fid = fopen(['Association' num2str(i) '.txt'],'w');
    for j=1:length(origRecon2.rxns)
        interactingGenes = origRecon2Genes( origRecon2.rxnGeneMat(j,:)~=0 );
        if length(interactingGenes) > 0
            fprintf(fid,'%s\t',origRecon2.rxns{j});
        else
            fprintf(fid,'%s\n',origRecon2.rxns{j});
        end
        for k=1:length(interactingGenes)-1
            fprintf(fid,'%s\t',interactingGenes{k});
        end
        if length(interactingGenes) > 0
            fprintf(fid,'%s\n',interactingGenes{end});
        end
    end
    fclose(fid);
end