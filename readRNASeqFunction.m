function readRNASeqFunction(rawData, rawMetadata)

load('readRNASeqDataMitoCarta.mat','pvals1corr','pvals2corr','pvals3corr','pvals4corr','pvals5corr','pvals6corr', ...
    'mouseIDToLineNum','mouseIDToQLength','mouseIDToTissue','mouseIDToSex', ...
    'pvals1t','pvals2t','pvals3t','pvals4t','pvals5t','pvals6t');

testLysosomal = 0;testMitoMiner = 0;testMitoCarta = 1;
testT = 0; testCorr = 1; 
testAll = 0;            %never use with testCorr, will output too many images
                        %doesn't matter though if we run lysosomal or
                        %mitochondrial

if(testLysosomal)
    inputFile = fopen('Ensembl Mouse List.csv');
end
if(testMitoMiner)
    inputFile = fopen('MitoMiner Mouse List.csv');
end
if(testMitoCarta)
    inputFile = fopen('MitoCarta Mouse List.csv');
end

lookEnsemblIDs = {};
line = fgetl(inputFile);
line = fgetl(inputFile); % skip header of inputFile
while line ~= -1
    words = strsplit(line,',');
    if(testLysosomal)
        if(length(words)==2) %check if actual ENSMUSG, might be empty
            lookEnsemblIDs{end+1} = words{2};
        end
    end
    if(testMitoMiner)
        temp = words{16};
        lookEnsemblIDs{end+1} = temp(2:end-1);
    end
    if(testMitoCarta)
        temp = words{11};
        lookEnsemblIDs{end+1} = temp(2:end-1);
    end
    line = fgetl(inputFile);
end
fclose(inputFile);

lookEnsemblIDToLineNum = containers.Map;
for i=1:size(rawData,1)
    if(any(strcmp(rawData{i,1}, lookEnsemblIDs)))
        lookEnsemblIDToLineNum(rawData{i,1}) = i;
    end
end

if(testT)
    pvals1 = pvals1t;
    pvals2 = pvals2t;
    pvals3 = pvals3t;
    pvals4 = pvals4t;
    pvals5 = pvals5t;
    pvals6 = pvals6t;
end
if(testCorr)
    pvals1 = pvals1corr;
    pvals2 = pvals2corr;
    pvals3 = pvals3corr;
    pvals4 = pvals4corr;
    pvals5 = pvals5corr;
    pvals6 = pvals6corr;
end

pvals1(isnan(pvals1)) = 1;
[sortPvals1 idxs1] = sort(pvals1);
[sortPvals2 idxs2] = sort(pvals2);
[sortPvals3 idxs3] = sort(pvals3);
[sortPvals4 idxs4] = sort(pvals4);
[sortPvals5 idxs5] = sort(pvals5);
[sortPvals6 idxs6] = sort(pvals6);

allEnsemblIDs = {};
for i=2:size(rawData,1)
    allEnsemblIDs{end+1} = rawData{i,1};
end
cutoff1 = FDRCutoff(sortPvals1)
cutoff2 = FDRCutoff(sortPvals2)
cutoff3 = FDRCutoff(sortPvals3)
cutoff4 = FDRCutoff(sortPvals4)
cutoff5 = FDRCutoff(sortPvals5)
cutoff6 = FDRCutoff(sortPvals6)
sigEnsemblIDs1 = allEnsemblIDs(idxs1(1:cutoff1));
sigEnsemblIDs2 = allEnsemblIDs(idxs2(1:cutoff2));
sigEnsemblIDs3 = allEnsemblIDs(idxs3(1:cutoff3));
sigEnsemblIDs4 = allEnsemblIDs(idxs4(1:cutoff4));
sigEnsemblIDs5 = allEnsemblIDs(idxs5(1:cutoff5));
sigEnsemblIDs6 = allEnsemblIDs(idxs6(1:cutoff6));

lookEnsemblIDs = keys(lookEnsemblIDToLineNum);

lookSigEnsemblIDs1 = intersect(lookEnsemblIDs, sigEnsemblIDs1);
lookSigEnsemblIDs2 = intersect(lookEnsemblIDs, sigEnsemblIDs2);
lookSigEnsemblIDs3 = intersect(lookEnsemblIDs, sigEnsemblIDs3);
lookSigEnsemblIDs4 = intersect(lookEnsemblIDs, sigEnsemblIDs4);
lookSigEnsemblIDs5 = intersect(lookEnsemblIDs, sigEnsemblIDs5);
lookSigEnsemblIDs6 = intersect(lookEnsemblIDs, sigEnsemblIDs6);

% temp={};
% for i=1:length(lookSigEnsemblIDs3)
%     insExpressionVals31 = [];
%     insExpressionVals32 = [];
%     for j=2:size(rawData,2)
%         if(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==20)
%             insExpressionVals31(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs3{i}),j };
%         elseif(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==175)
%             insExpressionVals32(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs3{i}),j };
%         end
%     end
%     if(median(insExpressionVals32)>=2*median(insExpressionVals31))
%         temp{end+1} = lookSigEnsemblIDs3{i};
%     end
% end
% lookSigEnsemblIDs3 = temp;
% temp={};
% for i=1:length(lookSigEnsemblIDs4)
%     insExpressionVals41 = [];
%     insExpressionVals42 = [];
%     for j=2:size(rawData,2)
%         if(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==20)
%             insExpressionVals41(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs4{i}),j };
%         elseif(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==175)
%             insExpressionVals42(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs4{i}),j };
%         end
%     end
%     if(median(insExpressionVals42)>=2*median(insExpressionVals41))
%         temp{end+1} = lookSigEnsemblIDs4{i};
%     end
% end
% lookSigEnsemblIDs4 = temp;
% temp={};
% for i=1:length(lookSigEnsemblIDs6)
%     insExpressionVals61 = [];
%     insExpressionVals62 = [];
%     for j=2:size(rawData,2)
%         if(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==20)
%             insExpressionVals61(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs6{i}),j };
%         elseif(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && ...
%             strcmp(mouseIDToSex(rawData{1,j}),'F') && mouseIDToQLength(rawData{1,j})==175)
%             insExpressionVals62(end+1) = rawData{ ensemblIDToLineNum(lookSigEnsemblIDs6{i}),j };
%         end
%     end
%     if(median(insExpressionVals62)>=2*median(insExpressionVals61))
%         temp{end+1} = lookSigEnsemblIDs6{i};
%     end
% end
% lookSigEnsemblIDs6 = temp;
%
fid1 = fopen('lookSigEnsemblIDs1.txt','w');
for i=1:length(lookSigEnsemblIDs1)
    fprintf(fid1,'%s\n',lookSigEnsemblIDs1{i});
end
fclose(fid1);
fid2 = fopen('lookSigEnsemblIDs2.txt','w');
for i=1:length(lookSigEnsemblIDs2)
    fprintf(fid2,'%s\n',lookSigEnsemblIDs2{i});
end
fclose(fid2);
fid3 = fopen('lookSigEnsemblIDs3.txt','w');
for i=1:length(lookSigEnsemblIDs3)
    fprintf(fid3,'%s\n',lookSigEnsemblIDs3{i});
end
fclose(fid3);
fid4 = fopen('lookSigEnsemblIDs4.txt','w');
for i=1:length(lookSigEnsemblIDs4)
    fprintf(fid4,'%s\n',lookSigEnsemblIDs4{i});
end
fclose(fid4);
fid5 = fopen('lookSigEnsemblIDs5.txt','w');
for i=1:length(lookSigEnsemblIDs5)
    fprintf(fid5,'%s\n',lookSigEnsemblIDs5{i});
end
fclose(fid5);
fid6 = fopen('lookSigEnsemblIDs6.txt','w');
for i=1:length(lookSigEnsemblIDs6)
    fprintf(fid6,'%s\n',lookSigEnsemblIDs6{i});
end
fclose(fid6);
fid1 = fopen('sigEnsemblIDs1.txt','w');
for i=1:length(sigEnsemblIDs1)
    fprintf(fid1,'%s\n',sigEnsemblIDs1{i});
end
fclose(fid1);
fid2 = fopen('sigEnsemblIDs2.txt','w');
for i=1:length(sigEnsemblIDs2)
    fprintf(fid2,'%s\n',sigEnsemblIDs2{i});
end
fclose(fid2);
fid3 = fopen('sigEnsemblIDs3.txt','w');
for i=1:length(sigEnsemblIDs3)
    fprintf(fid3,'%s\n',sigEnsemblIDs3{i});
end
fclose(fid3);
fid4 = fopen('sigEnsemblIDs4.txt','w');
for i=1:length(sigEnsemblIDs4)
    fprintf(fid4,'%s\n',sigEnsemblIDs4{i});
end
fclose(fid4);
fid5 = fopen('sigEnsemblIDs5.txt','w');
for i=1:length(sigEnsemblIDs5)
    fprintf(fid5,'%s\n',sigEnsemblIDs5{i});
end
fclose(fid5);
fid6 = fopen('sigEnsemblIDs6.txt','w');
for i=1:length(sigEnsemblIDs6)
    fprintf(fid6,'%s\n',sigEnsemblIDs6{i});
end
fclose(fid6);
fid1 = fopen('allEnsemblIDs.txt','w');
for i=1:length(allEnsemblIDs)
    fprintf(fid1,'%s\n',allEnsemblIDs{i});
end
fclose(fid1);
fid2 = fopen('lookEnsemblIDs.txt','w');
for i=1:length(lookEnsemblIDs)
    fprintf(fid2,'%s\n',lookEnsemblIDs{i});
end
fclose(fid2);

%length(lookEnsemblIDs)/37620*100
%sqrt(100*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620)
overallP1=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs1)-length(lookEnsemblIDs)/37620*cutoff1)/sqrt(cutoff1*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP2=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs2)-length(lookEnsemblIDs)/37620*cutoff2)/sqrt(cutoff2*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP3=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs3)-length(lookEnsemblIDs)/37620*cutoff3)/sqrt(cutoff3*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP4=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs4)-length(lookEnsemblIDs)/37620*cutoff4)/sqrt(cutoff4*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP5=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs5)-length(lookEnsemblIDs)/37620*cutoff5)/sqrt(cutoff5*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP6=1-cdf(makedist('Normal'),(length(lookSigEnsemblIDs6)-length(lookEnsemblIDs)/37620*cutoff6)/sqrt(cutoff6*length(lookEnsemblIDs)/37620*(37620-length(lookEnsemblIDs))/37620))
overallP1Hyper=1-hygecdf(length(lookSigEnsemblIDs1),37620,length(lookEnsemblIDs),length(sigEnsemblIDs1))
overallP2Hyper=1-hygecdf(length(lookSigEnsemblIDs2),37620,length(lookEnsemblIDs),length(sigEnsemblIDs2))
overallP3Hyper=1-hygecdf(length(lookSigEnsemblIDs3),37620,length(lookEnsemblIDs),length(sigEnsemblIDs3))
overallP4Hyper=1-hygecdf(length(lookSigEnsemblIDs4),37620,length(lookEnsemblIDs),length(sigEnsemblIDs4))
overallP5Hyper=1-hygecdf(length(lookSigEnsemblIDs5),37620,length(lookEnsemblIDs),length(sigEnsemblIDs5))
overallP6Hyper=1-hygecdf(length(lookSigEnsemblIDs6),37620,length(lookEnsemblIDs),length(sigEnsemblIDs6))
enrichmentMatrix = [length(sigEnsemblIDs1) length(lookSigEnsemblIDs1); ...
    length(sigEnsemblIDs2) length(lookSigEnsemblIDs2); ...
    length(sigEnsemblIDs3) length(lookSigEnsemblIDs3); ...
    length(sigEnsemblIDs4) length(lookSigEnsemblIDs4); ...
    length(sigEnsemblIDs5) length(lookSigEnsemblIDs5); ...
    length(sigEnsemblIDs6) length(lookSigEnsemblIDs6)];

if(testAll)
    lookSigEnsemblIDs1 = sigEnsemblIDs1;
    lookSigEnsemblIDs2 = sigEnsemblIDs2;
    lookSigEnsemblIDs3 = sigEnsemblIDs3;
    lookSigEnsemblIDs4 = sigEnsemblIDs4;
    lookSigEnsemblIDs5 = sigEnsemblIDs5;
    lookSigEnsemblIDs6 = sigEnsemblIDs6;
end

lookSigEnsemblVals1 = [];
lookSigEnsemblQLengths1 = [];
for i=1:length(lookSigEnsemblIDs1)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'Liver') && strcmp(mouseIDToSex(rawData{1,j}),'F'))
            %find(strcmp(lookSigEnsemblIDs2{i},allEnsemblIDs))+1
            lookSigEnsemblVals1(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs1{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths1(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end
lookSigEnsemblVals2 = [];
lookSigEnsemblQLengths2 = [];
for i=1:length(lookSigEnsemblIDs2)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'Liver') && strcmp(mouseIDToSex(rawData{1,j}),'F'))
            %find(strcmp(lookSigEnsemblIDs2{i},allEnsemblIDs))+1
            lookSigEnsemblVals2(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs2{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths2(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end
lookSigEnsemblVals3 = [];
lookSigEnsemblQLengths3 = [];
for i=1:length(lookSigEnsemblIDs3)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && strcmp(mouseIDToSex(rawData{1,j}),'F'))
            lookSigEnsemblVals3(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs3{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths3(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end
lookSigEnsemblVals4 = [];
lookSigEnsemblQLengths4 = [];
for i=1:length(lookSigEnsemblIDs4)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'cortex') && strcmp(mouseIDToSex(rawData{1,j}),'M'))
            lookSigEnsemblVals4(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs4{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths4(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end
lookSigEnsemblVals5 = [];
lookSigEnsemblQLengths5 = [];
for i=1:length(lookSigEnsemblIDs5)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'striatum') && strcmp(mouseIDToSex(rawData{1,j}),'M'))
            lookSigEnsemblVals5(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs5{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths5(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end
lookSigEnsemblVals6 = [];
lookSigEnsemblQLengths6 = [];
for i=1:length(lookSigEnsemblIDs6)
    count=1;
    for j=2:size(rawData,2)
        if(strcmp(mouseIDToTissue(rawData{1,j}),'Liver') && strcmp(mouseIDToSex(rawData{1,j}),'F'))
            %find(strcmp(lookSigEnsemblIDs2{i},allEnsemblIDs))+1
            lookSigEnsemblVals6(i,count) = rawData{ find(strcmp(lookSigEnsemblIDs6{i},allEnsemblIDs))+1,j };
            lookSigEnsemblQLengths6(count) = mouseIDToQLength(rawData{1,j});
            count=count+1;
        end
    end
end

% lookSigEnsemblQLengths = [];
% for i=2:size(rawData,2)
%     lookSigEnsemblQLengths(end+1) = mouseIDToQLength(rawData{1,i});
% end
% 
% lookSigEnsemblQLengths = transformQVals(lookSigEnsemblQLengths);
lookSigEnsemblQLengths1 = transformQVals(lookSigEnsemblQLengths1);
lookSigEnsemblQLengths2 = transformQVals(lookSigEnsemblQLengths2);
lookSigEnsemblQLengths3 = transformQVals(lookSigEnsemblQLengths3);
lookSigEnsemblQLengths4 = transformQVals(lookSigEnsemblQLengths4);
lookSigEnsemblQLengths5 = transformQVals(lookSigEnsemblQLengths5);
lookSigEnsemblQLengths6 = transformQVals(lookSigEnsemblQLengths6);

if(testLysosomal)
    if(testT)
        outputFile = 'lysosomalVariablesT.mat';
    end
    if(testCorr)
        outputFile = 'lysosomalVariablesCorr.mat';
    end
end
if(testMitoMiner)
    if(testT)
        outputFile = 'MitoMinerVariablesT.mat';
    end
    if(testCorr)
        outputFile = 'MitoMinerVariablesCorr.mat';
    end
end
if(testMitoCarta)
    if(testT)
        outputFile = 'MitoCartaVariablesT.mat';
    end
    if(testCorr)
        outputFile = 'MitoCartaVariablesCorr.mat';
    end
end
if(testAll)
    outputFile = 'AllVariablesT.mat';
end

save(outputFile, ...
    'lookSigEnsemblQLengths1','lookSigEnsemblQLengths5','lookSigEnsemblQLengths3','lookSigEnsemblQLengths4','lookSigEnsemblQLengths6','lookSigEnsemblQLengths2', ...
    'lookSigEnsemblVals1','lookSigEnsemblVals5','lookSigEnsemblVals3','lookSigEnsemblVals4','lookSigEnsemblVals6','lookSigEnsemblVals2', ...
    'lookSigEnsemblIDs1','lookSigEnsemblIDs5','lookSigEnsemblIDs2','lookSigEnsemblIDs3','lookSigEnsemblIDs4','lookSigEnsemblIDs6', ...
    'overallP1','overallP2','overallP3','overallP4','overallP5','overallP6','enrichmentMatrix');

    function cutoff = FDRCutoff(sortPvals)
        cutoff = 1;
        for i=1:length(sortPvals)
            if(sortPvals(i)<=i*.05/length(sortPvals))
                cutoff=cutoff+1;
            else
                break;
            end
        end
    end

    function newQVals = transformQVals(oldQVals)
        newQVals = zeros(size(oldQVals));
        for innerI=1:length(oldQVals)
            if(oldQVals(innerI)==20)
                newQVals(innerI)=1;
            elseif(oldQVals(innerI)==80)
                newQVals(innerI)=2;
            elseif(oldQVals(innerI)==92)
                newQVals(innerI)=3;
            elseif(oldQVals(innerI)==111)
                newQVals(innerI)=4;
            elseif(oldQVals(innerI)==140)
                newQVals(innerI)=5;
            elseif(oldQVals(innerI)==175)
                newQVals(innerI)=6;
            end
        end
    end
end