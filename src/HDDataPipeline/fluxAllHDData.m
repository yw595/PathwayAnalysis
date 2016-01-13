<<<<<<< HEAD
config;
v2struct();
disp('running fluxAllHDData');
outputDir1 = [outputDir filesep 'fluxAllHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

outputPrefixes = {'mCADRE'};%{'Normal', 'iMAT', 'GIMME', ...
%    'iMATMachado', 'GIMMEMachado','EFlux','GXFBA'};

fluxResults = {};
for i=1:length(outputPrefixes)
    for j=1:size(allExpressionData,2)
        for k=1:2
            try
                initCobraToolbox;
                
                FBAModel = changeObjective(origRecon2,'biomass_reaction');
                if k==1
                    constrainedModel = constrainMediumExc(initializeRecon2(FBAModel));
                else
                    constrainedModel = constrainMediumExc2(initializeRecon2(FBAModel));
                end
                FBASoln = optimizeCbModel(constrainedModel,1);
                v_fba = FBASoln.x;
                
                if strcmp(outputPrefixes{i},'RELATCH')
                    tempIDs = {};
                    for z =1:length(expressionIDsMachado)
                        tempIDs{z} = num2str(expressionIDsMachado(z));
                    end
                    RELATCHFluxes = call_RELATCH(constrainMediumExc(initializeRecon2(origRecon2)), ...
                        tempIDs, expressionDataMachado, {}, []);
                elseif strcmp(outputPrefixes{i},'MADE')
                    [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
                        readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
                    [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
                    intersectDataA = expressionDataMachado(intersectIdxsA);
                    intersectDataB = expressionDataMachadoRef(intersectIdxsB);
                    bounds_ref = struct();
                    bounds_ref.lb = origRecon2.lb;
                    bounds_ref.ub = origRecon2.ub;
                    MADEFluxes = call_MADE(constrainMediumExc(initializeRecon2(origRecon2)), ...
                        num2str(intersectIDs), intersectDataA, intersectDataB, 0.9, [], bounds_ref);
                elseif strcmp(outputPrefixes{i},'EFlux')
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
                else
                    runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
                        [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);
                    
                    modelToRun = addBiomass(origRecon2,modelToRun);
                    modelToRun = changeObjective(modelToRun,'biomass_reaction');
                    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(modelToRun)),1);
                    v_fba = FBASoln.x;
                end
            catch
                disp(['Error at ' num2str(i) ' ' num2str(j) ' ' num2str(k)]);
            end
        end
    end
end