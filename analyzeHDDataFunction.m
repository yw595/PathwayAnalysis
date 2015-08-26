function analyzeHDDataFunction()

analyzeAll=1;

timePoint=3;
if analyzeAll
    load('AllHDData.mat','allObservationIDs','allMouseIDs','allTissues','allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');
else
    load('HDData.mat','mouseIDs','tissues','QLengths','sexes','geneIDs','expressionDataHD');
end
if analyzeAll
    load('AllHDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');
    pValsCorr = squeeze(pValsCorr(timePoint,:,:));
else
    load('HDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');
end
testLysosomal = 0;testMitoMiner = 0;testMitoCarta = 1;
testT = 0; testCorr = 1;
%never use with testCorr, will output too many images
%doesn't matter though if we run lysosomal or mitochondrial
testAll = 0;            

if(testLysosomal)
    inputFile = fopen(['input' filesep 'Ensembl Mouse List.csv']);
    dataFields = textscan(inputFile,repmat('%s',1,2),'Delimiter',',','HeaderLines',1);
elseif(testMitoMiner)
    inputFile = fopen(['input' filesep 'MitoMiner Mouse Ensembl To Human Entrez.csv']);
    dataFields = textscan(inputFile,repmat('%s',1,16),'Delimiter',',','HeaderLines',1);
elseif(testMitoCarta)
    if analyzeAll
        inputFile = fopen(['input' filesep 'MitoCarta Mouse Ensembl To Mouse Entrez.csv']);
    else
        inputFile = fopen(['input' filesep 'MitoCarta Mouse Ensembl To Human Entrez.csv']);
    end
    dataFields = textscan(inputFile,repmat('%s',1,2),'Delimiter',',','HeaderLines',1);
end
fclose(inputFile);

lookEnsemblIDs = dataFields{1};
lookEntrezIDs = dataFields{2};
ensemblToEntrez = containers.Map(lookEnsemblIDs,lookEntrezIDs);

if(testT)
    pVals = pValsT;
end
if(testCorr)
    pVals = pValsCorr;
end

pVals(isnan(pVals)) = 1;
enrichmentMatrix = [];
for i=1:size(pVals,2)
    [sortPValsCol idxs] = sort(pVals(:,i));
    cutoffs(i) = FDRCutoff(sortPValsCol);

    if analyzeAll
        sigEntrezIDs{i} = allGeneIDs(idxs(1:cutoffs(i))); sortPValsCol = sortPValsCol(1:cutoffs(i));
        [lookSigEntrezIDs{i}, ~, intersectIdxs] = intersect(lookEntrezIDs, sigEntrezIDs{i});
        length(lookSigEntrezIDs{i})
        length(sigEntrezIDs{i})
        enrichmentMatrix(i,2) = length(lookSigEntrezIDs{i});
        enrichmentMatrix(i,1) = length(sigEntrezIDs{i});
    else
        sigEnsemblIDs{i} = geneIDs(idxs(1:cutoffs(i))); sortPValsCol = sortPValsCol(1:cutoffs(i));
        [lookSigEnsemblIDs{i}, ~, intersectIdxs] = intersect(lookEnsemblIDs, sigEnsemblIDs{i});
        length(lookSigEnsemblIDs{i})
        length(sigEntrezIDs{i})
        enrichmentMatrix(i,2) = length(lookSigEntrezIDs{i});
        enrichmentMatrix(i,1) = length(sigEntrezIDs{i});
    end
    
    sortPValsCol = sortPValsCol(intersectIdxs);
    
    if analyzeAll
    else
        lookSigEntrezIDs = {};
        for j=1:length(lookSigEnsemblIDs{i})
            lookSigEntrezIDs{j} = ensemblToEntrez(lookSigEnsemblIDs{i}{j});
        end
        length(lookSigEntrezIDs)
    end
    
    sortPValsColArray{i} = sortPValsCol;
    lookSigEntrezIDsArray{i} = lookSigEntrezIDs;
end

if testLysosomal
    subDir = 'Lysosomal';
elseif testMitoCarta
    subDir = 'MitoCarta';
elseif testMitoMiner
    subDir = 'MitoMiner';
end

if analyzeAll
    subDir = ['analyzeAll' subDir];
end

if testT
    subDir = [subDir 'T'];
elseif testCorr
    subDir = [subDir 'Corr'];
end
if ~exist(['input' filesep subDir], 'dir')
    system(['mkdir input' filesep subDir]);
end

for i=1:size(pVals,2)
    if analyzeAll
        fid = fopen(['input' filesep subDir filesep 'lookSigEntrezIDs' num2str(i) '.txt'],'w');
        ithLookSigEntrezIDs = lookSigEntrezIDs{i};
        for j=1:length(ithLookSigEntrezIDs)
            fprintf(fid,'%s\n',ithLookSigEntrezIDs{j});
        end
        fclose(fid);
    else
        fid = fopen(['input' filesep subDir filesep 'lookSigEnsemblIDs' num2str(i) '.txt'],'w');
        ithLookSigEnsemblIDs = lookSigEnsemblIDs{i};
        for j=1:length(ithLookSigEnsemblIDs)
            fprintf(fid,'%s\n',ithLookSigEnsemblIDs{j});
        end
        fclose(fid);
    end
end

for i=1:size(pVals,2)
    if analyzeAll
        fid = fopen(['input' filesep subDir filesep 'sigEntrezIDs' num2str(i) '.txt'],'w');
        ithSigEntrezIDs = sigEntrezIDs{i};
        for j=1:length(ithSigEntrezIDs)
            fprintf(fid,'%s\n',ithSigEntrezIDs{j});
        end
        fclose(fid);
    else
        fid = fopen(['input' filesep subDir filesep 'sigEnsemblIDs' num2str(i) '.txt'],'w');
        ithSigEnsemblIDs = sigEnsemblIDs{i};
        for j=1:length(ithSigEnsemblIDs)
            fprintf(fid,'%s\n',ithSigEnsemblIDs{j});
        end
        fclose(fid);
    end
end

for i=1:size(pVals,2)
    if analyzeAll
        overallPArray(i) = 1-cdf(makedist('Normal'),(length(lookSigEntrezIDs{i})-length(lookEntrezIDs)/23000*cutoffs(i)) ...
            /sqrt(cutoffs(i)*length(lookEntrezIDs)/23000*(23000-length(lookEntrezIDs))/23000));
        overallPHyperArray(i) = 1-hygecdf(length(lookSigEntrezIDs{i}),23000,length(lookEntrezIDs),length(sigEntrezIDs{i}));
    else
        overallPArray(i) = 1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs{i})-length(lookEnsemblIDs)/37620*cutoffs(i)) ...
            /sqrt(cutoffs(i)*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620));
        overallPHyperArray(i) = 1-hygecdf(length(lookSigEnsemblIDs{i}),37620,length(lookEnsemblIDs),length(sigEnsemblIDs{i}));
    end
end
overallPArray
overallPHyperArray

if(testAll)
    if analyzeAll
        for i=1:size(pVals,2)
            lookSigEntrezIDs{i} = sigEntrezIDs{i};
        end
    else
        for i=1:size(pVals,2)
            lookSigEnsemblIDs{i} = sigEnsemblIDs{i};
        end
    end
end

count = 0;
for k=1:length(uniqueTissues)
    if analyzeAll
        count = count+1;
        ithLookSigEntrezIDs = lookSigEntrezIDs{count};
        if testCorr
            mask = strcmp(allTissues, uniqueTissues{k}) ...
            & allMonths==allMonths(timePoint) & strcmp(allSeqTypes,'mRNA');
        else
            mask = strcmp(allTissues, uniqueTissues{k}) & (allQLengths==20 | allQLengths==175) ...
            & (allMonths==2 | allMonths==10) & strcmp(allSeqTypes,'mRNA');
        end
        [~, intersectIdxs, ~] = intersect(ithLookSigEntrezIDs,allGeneIDs);
        lookSigEntrezVals{count} = allExpressionData( intersectIdxs,mask );
        lookSigEntrezQLengths{count} = allQLengths(mask);
    else
        for l=1:length(uniqueSexes)
            count = count+1;
            ithLookSigEnsemblIDs = lookSigEnsemblIDs{count};
            mask = strcmp(tissues, uniqueTissues{k}) & strcmp(sexes, uniqueSexes{l});
            [~, intersectIdxs, ~] = intersect(ithLookSigEnsemblIDs,geneIDs);
            lookSigEnsemblVals{count} = expressionDataHD( intersectIdxs,mask );
            lookSigEnsemblQLengths{count} = cellfun(@(x) str2num(x(2:end)), QLengths(mask));
        end
    end
end

if analyzeAll
    for i=1:length(lookSigEntrezQLengths)
        lookSigEntrezQLengths{i} = transformQVals(lookSigEntrezQLengths{i});
    end
else
    for i=1:length(lookSigEnsemblQLengths)
        lookSigEnsemblQLengths{i} = transformQVals(lookSigEnsemblQLengths{i});
    end
end

outputFile = ['output' filesep 'MATFiles'];
if testLysosomal
    outputFile = [outputFile 'lysosomalVariables'];
elseif testMitoMiner
    outputFile = [outputFile 'MitoMinerVariables'];
elseif testMitoCarta
    outputFile = [outputFile 'MitoCartaVariables'];
elseif testAll
    outputFile = [outputFile 'AllVariables'];
end

if testT
    outputFile = [outputFile 'T.mat'];
elseif testCorr
    outputFile = [outputFile 'Corr.mat'];
end

if analyzeAll
    save(outputFile, 'lookSigEntrezQLengths', 'lookSigEntrezVals', 'lookSigEntrezIDs', 'overallPArray', 'overallPHyperArray',...
        'sortPValsColArray','enrichmentMatrix');
else
    save(outputFile, 'lookSigEnsemblQLengths', 'lookSigEnsemblVals', 'lookSigEnsemblIDs', 'overallPArray', 'overallPHyperArray',...
        'sortPValsColArray','lookSigEntrezIDsArray','enrichmentMatrix');
end

    function cutoff = FDRCutoff(sortPvals)
        cutoff = 1;
        for index1=1:length(sortPvals)
            if(sortPvals(index1)<=index1*.05/length(sortPvals))
                cutoff=cutoff+1;
            else
                break;
            end
        end
    end

    function newQVals = transformQVals(oldQVals)
        newQVals = zeros(size(oldQVals));
        if analyzeAll
            map = [20 1;80 0;92 0;111 0;140 0;175 2];
        else
            map = [20 1;80 2;92 3;111 4;140 5;175 6];
        end
        for index2=1:size(map,1)
            newQVals(oldQVals==map(index2,1))=map(index2,2);
        end
    end
end