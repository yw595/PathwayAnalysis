function analyzeHDDataFunction()

load('HDData.mat','mouseIDs','tissues','QLengths','sexes','geneIDs','expressionDataHD');
load('HDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');
testLysosomal = 0;testMitoMiner = 0;testMitoCarta = 1;
testT = 0; testCorr = 1; 
testAll = 0;            %never use with testCorr, will output too many images
                        %doesn't matter though if we run lysosomal or mitochondrial

if(testLysosomal)
    inputFile = fopen(['input' filesep 'Ensembl Mouse List.csv']);
    dataFields = textscan(inputFile,repmat('%s',1,2),'Delimiter',',','HeaderLines',1);
elseif(testMitoMiner)
    inputFile = fopen(['input' filesep 'MitoMiner Mouse Ensembl To Human Entrez.csv']);
    dataFields = textscan(inputFile,repmat('%s',1,16),'Delimiter',',','HeaderLines',1);
elseif(testMitoCarta)
    inputFile = fopen(['input' filesep 'MitoCarta Mouse Ensembl To Human Entrez.csv']);
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
for i=1:size(pVals,2)
    [sortPValsCol idxs] = sort(pVals(:,i));
    cutoffs(i) = FDRCutoff(sortPValsCol);
    sigEnsemblIDs{i} = geneIDs(idxs(1:cutoffs(i))); sortPValsCol = sortPValsCol(1:cutoffs(i));
    [lookSigEnsemblIDs{i}, ~, intersectIdxs] = intersect(lookEnsemblIDs, sigEnsemblIDs{i});
    length(lookSigEnsemblIDs{i})
    length(intersectIdxs)
    sortPValsCol = sortPValsCol(intersectIdxs);
    lookSigEntrezIDs = {};
    for j=1:length(lookSigEnsemblIDs{i})
        lookSigEntrezIDs{j} = ensemblToEntrez(lookSigEnsemblIDs{i}{j});
    end
    length(lookSigEntrezIDs)
    sortPValsColArray{i} = sortPValsCol;
    lookSigEntrezIDsArray{i} = lookSigEntrezIDs;
end

for i=1:size(pVals,2)
    fid = fopen(['output' filesep 'lookSigEnsemblIDs' num2str(i) '.txt'],'w');
    ithLookSigEnsemblIDs = lookSigEnsemblIDs{i};
    for j=1:length(ithLookSigEnsemblIDs)
        fprintf(fid,'%s\n',ithLookSigEnsemblIDs{j});
    end
    fclose(fid);
end

for i=1:size(pVals,2)
    fid = fopen(['output' filesep 'sigEnsemblIDs' num2str(i) '.txt'],'w');
    ithSigEnsemblIDs = sigEnsemblIDs{i};
    for j=1:length(ithSigEnsemblIDs)
        fprintf(fid,'%s\n',ithSigEnsemblIDs{j});
    end
    fclose(fid);
end

for i=1:size(pVals,2)
    overallPArray(i) = 1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs{i})-length(lookEnsemblIDs)/37620*cutoffs(i)) ...
        /sqrt(cutoffs(i)*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620));
    overallPHyperArray(i) = 1-hygecdf(length(lookSigEnsemblIDs{i}),37620,length(lookEnsemblIDs),length(sigEnsemblIDs{i}));
end

if(testAll)
    for i=1:size(pVals,2)
        lookSigEnsemblIDs{i} = sigEnsemblIDs{i};
    end
end

count = 0;
for k=1:length(uniqueTissues)
    for l=1:length(uniqueSexes)
        count = count+1;
        ithLookSigEnsemblIDs = lookSigEnsemblIDs{count};
        mask = strcmp(tissues, uniqueTissues{k}) & strcmp(sexes, uniqueSexes{l});
        [~, intersectIdxs, ~] = intersect(ithLookSigEnsemblIDs,geneIDs);
        lookSigEnsemblVals{count} = expressionDataHD( intersectIdxs,mask );
        lookSigEnsemblQLengths{count} = cellfun(@(x) str2num(x(2:end)), QLengths(intersectIdxs));
    end
end

for i=1:length(lookSigEnsemblQLengths)
    lookSigEnsemblQLengths{i} = transformQVals(lookSigEnsemblQLengths{i});
end

outputFile = ['output' filesep];
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

save(outputFile, 'lookSigEnsemblQLengths', 'lookSigEnsemblVals', 'lookSigEnsemblIDs', 'overallPArray', ...
    'sortPValsColArray','lookSigEntrezIDsArray');

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
        map = [20 1;80 2;92 3;111 4;140 5;175 6];
        for index2=1:size(map,1)
            newQVals(oldQVals==map(index2,1))=map(index2,2);
        end
    end
end