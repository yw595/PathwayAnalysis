configPicard;
v2struct();
disp('running processPicardData');
outputDir1 = [outputDir filesep 'processPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

labelsT = {'heteroplasmy'};
labelsCorr = {'heteroplasmy'};
mask1 = strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{1});
mask2 = strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{end-1});
[~, pValsT] = ttest2(PicardExpressionData(:,mask1)', PicardExpressionData(:,mask2)');
mask1 = ~strcmp(PicardHeteroplasmies,'Rho0');
[rhos pValsCorr] = corr(cellfun(@(x) str2num(x), PicardHeteroplasmies(mask1)), PicardExpressionData(:,mask1)', 'type', 'Spearman');

save([outputDir1 filesep 'PicardDataProcessed.mat'],'pValsCorr','pValsT','rhos','labelsT','labelsCorr');
