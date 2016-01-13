configPicard;
v2struct();
disp('running RMAPicardData');
outputDir1 = [outputDir filesep 'RMAPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if ~exist([outputDir1 filesep 'anovaPVals.mat'],'file')
    anovaPVals = arrayfun(@(n) anova1(PicardExpressionData(n,:), PicardHeteroplasmies,'off'),1:size(PicardExpressionData,1));
    save([outputDir1 filesep 'anovaPVals.mat',anovaPVals);
else
    load([outputDir1 filesep 'anovaPVals.mat']);
end

[ensemblPicardGeneIDs pValsT pValsCorr] = mapGeneIDs([inputDir filesep 'humanEnsemblToHumanHGNC.csv'],PicardGeneIDs,1, 'pVals1', pValsT, 'pVals2', pValsCorr);
%[entrezPicardGeneIDs pValsT pValsCorr] = mapGeneIDs([inputDir filesep 'humanEnsemblToHumanEntrez.csv'],ensemblPicardGeneIDs,0, 'pVals1', pValsT, 'pVals2', pValsCorr);
[entrezPicardGeneIDs pValsAnova] = mapGeneIDs([inputDir filesep 'humanEnsemblToHumanEntrez.csv'],ensemblPicardGeneIDs,0, 'pVals1', anovaPVals);

%[sortedNormScoresT sortedMetsT sortedMetNamesT] = RMAPicardDataFunction(entrezPicardGeneIDs, pValsT, origRecon2, 1, 0);
%[sortedNormScoresCorr sortedMetsCorr sortedMetNamesCorr] = RMAPicardDataFunction(entrezPicardGeneIDs, pValsCorr, origRecon2, 1, 0);
[sortedNormScoresAnova sortedMetsAnova sortedMetNamesAnova] = RMAPicardDataFunction(entrezPicardGeneIDs, pValsAnova, origRecon2, 0, 0);

%refNormScoresT = sortedNormScoresT; refNormScoresCorr = sortedNormScoresCorr;
%refMetsT = sortedMetsT; refMetsCorr = sortedMetsCorr;
%refMetNamesT = sortedMetNamesT; refMetNamesCorr = sortedMetNamesCorr;
%refPValsT = pValsT; refPValsCorr = pValsCorr;
refNormScoresAnova = sortedNormScoresAnova; refMetsAnova = sortedMetsAnova; refMetNamesAnova = sortedMetNamesAnova; refPValsAnova = pValsAnova;

for j=1:1%2
    for i=1:10
        disp([num2str(j) ' ' num2str(i)]);

        if j==1
            sigCutoff = 0;
            permuteCutoff = 1;
        else
            sigCutoff = 1;
            permuteCutoff = 2;
        end
        
	[sortedNormScoresAnova sortedMetsAnova sortedMetNamesAnova] = RMAPicardDataFunction(entrezPicardGeneIDs, refPValsAnova, origRecon2, sigCutoff, permuteCutoff);
        [~, intersectIdxs] = ismember(sortedMetsAnova, refMetsAnova);
        permuteScoresMatrixAnova(:,i,j) = sortedNormScoresAnova(intersectIdxs);

        %[sortedNormScoresT sortedMetsT sortedMetNamesT] = RMAPicardDataFunction(entrezPicardGeneIDs, refPValsT, origRecon2, sigCutoff, permuteCutoff);
        %[~, intersectIdxs] = ismember(sortedMetsT, refMetsT);
        %permuteScoresMatrixT(:,i,j) = sortedNormScoresT(intersectIdxs);
        
        %[sortedNormScoresCorr sortedMetsCorr sortedMetNamesCorr] =RMAPicardDataFunction(entrezPicardGeneIDs, refPValsCorr,origRecon2, sigCutoff, permuteCutoff);
        %[~, intersectIdxs] = ismember(sortedMetsCorr, refMetsCorr);
        %permuteScoresMatrixCorr(:,i,j) = sortedNormScoresCorr(intersectIdxs);
    end
end

tempNames1 = {'refNormScores','refMets','refMetNames','refPVals', 'permuteScoresMatrix'};
tempNames2 = {'T','Corr','Anova'};
first = 1;
for i=1:length(tempNames1)
    for j=2:length(tempNames2)
        varName = [tempNames1{i} tempNames2{j}];
        if exist(varName,'var')
            if first
                save([outputDir1 filesep 'RMAPicardDataPermute.mat',varName);
                first=0;
            else
                save([outputDir1 filesep 'RMAPicardDataPermute.mat',varName,'append');
            end
        end
    end
end


%corr(cellfun(@(x) length(x),
%sortedSubsystemCorr),sortedNormScoresCorr,'type','Pearson')
%sortedSubsystemCorr = cellfun(@(x)
%origRecon2.subSystems(origRecon2.S(strcmp(x,origRecon2.mets),:)~=0),
%sortedMetsCorr,'UniformOutput',0)
