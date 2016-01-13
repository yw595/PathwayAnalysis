configPicard;
v2struct();
disp('running inputFALCONPicardData');
outputDir1 = [outputDir filesep 'inputFALCONPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

[ensemblPicardGeneIDs pValsT pValsCorr PicardExpressionData] = mapGeneIDs ...
    ([inputDir filesep 'humanEnsemblToHumanHGNC.csv'],PicardGeneIDs,1, 'pVals1', pValsT, 'pVals2' , pValsCorr, 'expressionData', PicardExpressionData);
[entrezPicardGeneIDs pValsT pValsCorr PicardExpressionData] = mapGeneIDs ...
    ([inputDir filesep 'humanEnsemblToHumanEntrez.csv'],ensemblPicardGeneIDs,0, 'pVals1', pValsT, 'pVals2' , pValsCorr, 'expressionData', PicardExpressionData);

for i=1:length(uniqueHeteroplasmies)
    disp(uniqueHeteroplasmies{i})
    relevantData = PicardExpressionData(:, strcmp(PicardHeteroplasmies,uniqueHeteroplasmies{i}));
    
    meanData = mean(relevantData,2);
    stdData = std(relevantData,0,2);
    
    outputFile = fopen([outputDir1 filesep 'FALCONData' num2str(uniqueHeteroplasmies{i}) '.txt'],'w');
    %fprintf(outputFile,'entrez Gene ID\tmean\tstd\n');
    for j=1:length(entrezPicardGeneIDs)
        fprintf(outputFile,'%s\t%f\t%f',entrezPicardGeneIDs{j}, meanData(j), stdData(j));
        if j~=length(entrezPicardGeneIDs)
            fprintf(outputFile,'\n');
        end
    end
    fclose(outputFile);
end