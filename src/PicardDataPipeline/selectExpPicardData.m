function selectExpPicardData(hugoGeneIDs, configStruct, anOutputDir)

v2struct(configStruct);
count=0;
xlabels = uniqueHeteroplasmies;

for i=1:length(hugoGeneIDs)
    disp(hugoGeneIDs{i});
    count=0;
    vals = [];
    xvals = [];
    for j=1:length(uniqueHeteroplasmies)
        count=count+1;
        disp([num2str(i) ' ' num2str(j)])
        [~, intersectIdxs1, intersectIdxs2] = intersect(PicardGeneIDs, hugoGeneIDs{i});
        mask = strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{j});
        tempVals = PicardExpressionData(intersectIdxs1,mask);
        vals = [vals(:); tempVals(:)];
        xvals = [xvals(:); count*ones(length(tempVals(:)),1)];
    end
    makeScatter(xvals, vals, hugoGeneIDs{i}, xlabels, 'Expression', anOutputDir);
end

end