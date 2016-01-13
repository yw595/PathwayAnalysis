configPicard;
v2struct();
disp('running fluxPicardData');
outputDir1 = [outputDir filesep 'fluxPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

[ensemblPicardGeneIDs pValsT pValsCorr PicardExpressionData] = mapGeneIDs ...
    ([inputDir filesep 'humanEnsemblToHumanHGNC.csv'],PicardGeneIDs,1, 'pVals1', pValsT, 'pVals2' , pValsCorr, 'expressionData', PicardExpressionData);
[entrezPicardGeneIDs pValsT pValsCorr PicardExpressionData] = mapGeneIDs ...
    ([inputDir filesep 'humanEnsemblToHumanEntrez.csv'],ensemblPicardGeneIDs,0, 'pVals1', pValsT, 'pVals2' , pValsCorr, 'expressionData', PicardExpressionData);

outputPrefixes = {'FALCON','EFlux','iMAT','GIMME'};%,'GXFBA'};%, 'iMAT', 'GIMME', ...
%    'iMATMachado', 'GIMMEMachado','EFlux','GXFBA'};

initCobraToolbox;

%load([topDir filesep 'fluxResults.mat']);
%fluxResults = {};
firstHeteroplasmies = zeros(length(PicardHeteroplasmies),1);
for j=1:length(uniqueHeteroplasmies)
    for k=1:length(PicardHeteroplasmies)
        if strcmp(PicardHeteroplasmies{k},uniqueHeteroplasmies{j})
            firstHeteroplasmies(k)=1;
            break;
        end
    end
end
% for i=2:length(outputPrefixes)
% %     if strcmp(outputPrefixes{i},'FALCON')
% %         currentHeteroplasmy = '';
% %     end
%     for k=1:2
% 
%         for j=1:size(PicardExpressionData,2)
%             %try
%                 
%                 FBAModel = changeObjective(origRecon2,'biomass_reaction');
%                 if k==1
%                     constrainedModel = constrainMediumExc(initializeRecon2(FBAModel));
%                 else
%                     constrainedModel = constrainMediumExc2(initializeRecon2(FBAModel));
%                 end
%                 FBASoln = optimizeCbModel(constrainedModel,1);
%                 v_fba = FBASoln.x;
%                 
% %                 if strcmp(outputPrefixes{i},'RELATCH')
% %                     tempIDs = {};
% %                     for z =1:length(expressionIDsMachado)
% %                         tempIDs{z} = num2str(expressionIDsMachado(z));
% %                     end
% %                     RELATCHFluxes = call_RELATCH(constrainMediumExc(initializeRecon2(origRecon2)), ...
% %                         tempIDs, expressionDataMachado, {}, []);
% %                 elseif strcmp(outputPrefixes{i},'MADE')
% %                     [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
% %                         readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
% %                     [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
% %                     intersectDataA = expressionDataMachado(intersectIdxsA);
% %                     intersectDataB = expressionDataMachadoRef(intersectIdxsB);
% %                     bounds_ref = struct();
% %                     bounds_ref.lb = origRecon2.lb;
% %                     bounds_ref.ub = origRecon2.ub;
% %                     MADEFluxes = call_MADE(constrainMediumExc(initializeRecon2(origRecon2)), ...
% %                         num2str(intersectIDs), intersectDataA, intersectDataB, 0.9, [], bounds_ref);
%                 if strcmp(outputPrefixes{i},'EFlux')
%                     fluxResults{i,j,k} = call_EFlux(constrainedModel, ...
%                         entrezPicardGeneIDs, PicardExpressionData(:,j), 'EX_glc(e)',1);
%                 elseif strcmp(outputPrefixes{i},'iMAT')
%                     fluxResults{i,j,k} = call_iMAT(constrainedModel, ...
%                         entrezPicardGeneIDs, PicardExpressionData(:,j), prctile(PicardExpressionData(:,j),25), prctile(PicardExpressionData(:,j),75),1);
%                 elseif strcmp(outputPrefixes{i},'GIMME')
%                     fluxResults{i,j,k} = call_GIMME(constrainedModel, ...
%                         entrezPicardGeneIDs, PicardExpressionData(:,j), .9, prctile(PicardExpressionData(:,j),25));
%                 elseif strcmp(outputPrefixes{i},'GXFBA')
%                     GXFBAWT = struct(); GXFBAWT.wt_sol.x = v_fba;
%                     GXFBAWT.wt_minf = origRecon2FVAMin; GXFBAWT.wt_maxf = origRecon2FVAMax;
%                     fluxResults{i,j,k} = call_GXFBA(constrainedModel, ...
%                         entrezPicardGeneIDs, PicardExpressionData(:,j), PicardExpressionData(:,1), GXFBAWT);
% %                 elseif strcmp(outputPrefixes{i},'FALCON')
% %                     if strcmp(PicardHeteroplasmies{j},currentHeteroplasmy)
% %                         fluxResults{i,j,k} = fluxResults{i,j-1,k};
% %                     else
% %                         if k==2
% %                             currentHeteroplasmy = PicardHeteroplasmies{j};
% %                         end
% %                         fluxResults{i,j,k} = runFALCONStrippedPicard(constrainedModel, ...
% %                             [outputDir filesep 'inputFALCONPicardData' filesep 'FALCONData' PicardHeteroplasmies{j} '.txt'], 1, outputDir1);
% %                     end
%                 elseif strcmp(outputPrefixes{i},'FALCON')
%                     if firstHeteroplasmies(j)==1
%                         fluxResults{i,j,k} = runFALCONStrippedPicard(constrainedModel, ...
%                             [outputDir filesep 'inputFALCONPicardData' filesep 'FALCONData' PicardHeteroplasmies{j} '.txt'], 1, outputDir1);
%                     end
%                 end
%             %catch
%                 %disp(['Error at ' num2str(i) ' ' num2str(j) ' ' num2str(k)]);
%             %end
%         end
%     end
% end

% FALCONIdx = find(strcmp(outputPrefixes,'FALCON'));
% for j=1:size(PicardExpressionData,2)
%     if firstHeteroplasmies(j)==0
%         prevIdx = j;
%         while firstHeteroplasmies(prevIdx)==0
%             prevIdx = prevIdx-1;
%         end
%         for k=1:2
%             fluxResults{FALCONIdx,j,k} = fluxResults{FALCONIdx,prevIdx,k};
%         end
%     end
% end

permuteFluxResults = {};

for m=1:2
    PicardExpressionData = PicardExpressionData(randsample(length(PicardExpressionData),length(PicardExpressionData)),:);
    
    fluxResults = {};
    for i=1:length(outputPrefixes)
        %     if strcmp(outputPrefixes{i},'FALCON')
        %         currentHeteroplasmy = '';
        %     end
        for k=1:2
            
            for j=1:size(PicardExpressionData,2)
                disp([num2str(m) ' ' num2str(i) ' ' num2str(k) ' ' num2str(j)]);
                try
                
                FBAModel = changeObjective(origRecon2,'biomass_reaction');
                if k==1
                    constrainedModel = constrainMediumExc(initializeRecon2(FBAModel));
                else
                    constrainedModel = constrainMediumExc2(initializeRecon2(FBAModel));
                end
                FBASoln = optimizeCbModel(constrainedModel,1);
                v_fba = FBASoln.x;
                
                %                 if strcmp(outputPrefixes{i},'RELATCH')
                %                     tempIDs = {};
                %                     for z =1:length(expressionIDsMachado)
                %                         tempIDs{z} = num2str(expressionIDsMachado(z));
                %                     end
                %                     RELATCHFluxes = call_RELATCH(constrainMediumExc(initializeRecon2(origRecon2)), ...
                %                         tempIDs, expressionDataMachado, {}, []);
                %                 elseif strcmp(outputPrefixes{i},'MADE')
                %                     [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
                %                         readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
                %                     [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
                %                     intersectDataA = expressionDataMachado(intersectIdxsA);
                %                     intersectDataB = expressionDataMachadoRef(intersectIdxsB);
                %                     bounds_ref = struct();
                %                     bounds_ref.lb = origRecon2.lb;
                %                     bounds_ref.ub = origRecon2.ub;
                %                     MADEFluxes = call_MADE(constrainMediumExc(initializeRecon2(origRecon2)), ...
                %                         num2str(intersectIDs), intersectDataA, intersectDataB, 0.9, [], bounds_ref);
                if strcmp(outputPrefixes{i},'EFlux')
                    fluxResults{i,j,k} = call_EFlux(constrainedModel, ...
                        entrezPicardGeneIDs, PicardExpressionData(:,j), 'EX_glc(e)',1);
                elseif strcmp(outputPrefixes{i},'iMAT')
                    fluxResults{i,j,k} = call_iMAT(constrainedModel, ...
                        entrezPicardGeneIDs, PicardExpressionData(:,j), prctile(PicardExpressionData(:,j),25), prctile(PicardExpressionData(:,j),75),1);
                elseif strcmp(outputPrefixes{i},'GIMME')
                    fluxResults{i,j,k} = call_GIMME(constrainedModel, ...
                        entrezPicardGeneIDs, PicardExpressionData(:,j), .9, prctile(PicardExpressionData(:,j),25));
                elseif strcmp(outputPrefixes{i},'GXFBA')
                    GXFBAWT = struct(); GXFBAWT.wt_sol.x = v_fba;
                    GXFBAWT.wt_minf = origRecon2FVAMin; GXFBAWT.wt_maxf = origRecon2FVAMax;
                    fluxResults{i,j,k} = call_GXFBA(constrainedModel, ...
                        entrezPicardGeneIDs, PicardExpressionData(:,j), PicardExpressionData(:,1), GXFBAWT);
                    %                 elseif strcmp(outputPrefixes{i},'FALCON')
                    %                     if strcmp(PicardHeteroplasmies{j},currentHeteroplasmy)
                    %                         fluxResults{i,j,k} = fluxResults{i,j-1,k};
                    %                     else
                    %                         if k==2
                    %                             currentHeteroplasmy = PicardHeteroplasmies{j};
                    %                         end
                    %                         fluxResults{i,j,k} = runFALCONStrippedPicard(constrainedModel, ...
                    %                             [outputDir filesep 'inputFALCONPicardData' filesep 'FALCONData' PicardHeteroplasmies{j} '.txt'], 1, outputDir1);
                    %                     end
                elseif strcmp(outputPrefixes{i},'FALCON')
                    if firstHeteroplasmies(j)==1
                        fluxResults{i,j,k} = runFALCONStrippedPicard(constrainedModel, ...
                            [outputDir filesep 'inputFALCONPicardData' filesep 'FALCONData' PicardHeteroplasmies{j} '.txt'], 1, outputDir1);
                    end
                end
                catch
                    disp(['Error at ' num2str(i) ' ' num2str(j) ' ' num2str(k)]);
                end
            end
        end
    end
    
    FALCONIdx = find(strcmp(outputPrefixes,'FALCON'));
    for j=1:size(PicardExpressionData,2)
        if firstHeteroplasmies(j)==0
            prevIdx = j;
            while firstHeteroplasmies(prevIdx)==0
                prevIdx = prevIdx-1;
            end
            for k=1:2
                fluxResults{FALCONIdx,j,k} = fluxResults{FALCONIdx,prevIdx,k};
            end
        end
    end
    
    permuteFluxResults{m} = fluxResults;
end
