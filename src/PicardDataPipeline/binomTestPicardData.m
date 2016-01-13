function pValsBinom = binomTestPicardData(hugoGeneIDs, configStruct)

v2struct(configStruct);

pValsBinom = zeros(length(uniqueHeteroplasmies),length(uniqueHeteroplasmies));
for i=1:length(uniqueHeteroplasmies)
    for j=1:length(uniqueHeteroplasmies)
        upDownReg = zeros(1,2);
        mask1 = strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{i});
        mask2 = strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{j});
        for k=1:length(hugoGeneIDs)
            [~, intersectIdxs1, intersectIdxs2] = intersect(PicardGeneIDs, hugoGeneIDs{k});
            baseVals = PicardExpressionData(intersectIdxs1,mask1);
            baseMedian = median(median(baseVals));
            diffVals = PicardExpressionData(intersectIdxs1,mask2);
            diffMedian = median(median(diffVals));
            if ttest2(baseVals,diffVals) < .05
                if diffMedian > baseMedian
                    upDownReg(1,1) = upDownReg(1,1) + 1;
                else
                    upDownReg(1,2) = upDownReg(1,2) + 1;
                end
            end
        end
        pValsBinom(i,j) = binocdf(upDownReg(1,1),sum(upDownReg(1,:),2),.5);
    end
end

end