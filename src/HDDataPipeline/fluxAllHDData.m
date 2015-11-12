outputPrefixes = {'mCADRE'};%{'Normal', 'iMAT', 'GIMME', ...
%    'iMATMachado', 'GIMMEMachado','EFlux','GXFBA'};

for i=1:length(outputPrefixes)
    initCobraToolbox;

    FBAModel = changeObjective(origRecon2,'biomass_reaction');
    FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(FBAModel)),1);
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
    elseif strcmp(outputPrefixes{j},'EFlux')
        EFluxes = call_EFlux(constrainMediumExc(initializeRecon2(origRecon2)), ...
            expressionIDsMachado, expressionDataMachado, 'EX_glc(e)',1);
    elseif strcmp(outputPrefixes{i},'GXFBA')
        [expressionIDsMachadoRef, expressionDataMachadoRef, ~] = ...
            readExpressionFile([inputDir filesep cellLinesArray{1} '.csv']);
        [intersectIDs, intersectIdxsA, intersectIdxsB] = intersect(expressionIDsMachado, expressionIDsMachadoRef);
        intersectDataA = expressionDataMachado(intersectIdxsA);
        intersectDataB = expressionDataMachadoRef(intersectIdxsB);
        
        GXFBAFluxes = call_GXFBA(constrainMediumExc(initializeRecon2(origRecon2)), ...
            num2str(intersectIDs), intersectDataA, intersectDataB, GXFBA);
    else
        runFALCONStripped(constrainMediumExc(initializeRecon2(modelToRun)), ...
            [inputDir filesep cellLinesArray{k} '.csv'], 1, outputDir);

        modelToRun = addBiomass(origRecon2,modelToRun);
        modelToRun = changeObjective(modelToRun,'biomass_reaction');
        FBASoln = optimizeCbModel(constrainMediumExc(initializeRecon2(modelToRun)),1);
        v_fba = FBASoln.x;
    end
end