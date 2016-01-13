config;
v2struct();
disp('running readCAndHHDData');
outputDir1 = [outputDir filesep 'readCAndHHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if isunix
    rawCAndHMatrix = readVaryFile([inputDir filesep 'cAndH' filesep 'GSE73508_series_matrix.txt'],'\t');
else
    [~, ~, rawCAndHMatrix] = xlsread([inputDir filesep 'cAndH' filesep 'GSE73508_series_matrix.xlsx']);
end

cAndHObservationIDs = cellfun(@(x) x(2:end-1), rawCAndHMatrix(50,2:end),'UniformOutput',0);
rawQLengths = rawCAndHMatrix(30,2:end);
for i=1:length(rawQLengths)
    if isempty(regexp(rawQLengths{i},'WT'))
        x = rawQLengths{i};
        cAndHQLengths{i} = [x(regexp(x,'(')+2:regexp(x,')')-1) 'Q'];
    else
        cAndHQLengths{i} = '0Q';
    end
end
cAndHMonths = cellfun(@(x) [x(regexp(x,':')+2:regexp(x,'mon')-2) 'Months'], rawCAndHMatrix(41,2:end), 'UniformOutput', 0);
%%%% EXTREMELY BAD ERROR, ROW DATA SWITCHED AT 501!!!!
cAndHTissues = cellfun(@(x) lower(x(regexp(x,'-')+2:end-1)), [rawCAndHMatrix(43,2:501) rawCAndHMatrix(44,502:end)],'UniformOutput',0);
for i=2:size(rawCAndHMatrix,2)
    if ~isempty(regexp(rawCAndHMatrix{30,i},'mRNA'));
        cAndHSeqTypes{i-1} = 'mRNA';
    else
        cAndHSeqTypes{i-1} = 'miRNA';
    end
end
cAndHSexes = cellfun(@(x) upper(x(regexp(x,':')+2:end)), rawCAndHMatrix(42,2:end),'UniformOutput',0);

cAndHExpressionData = NaN(39178,length(cAndHObservationIDs));

fileList = dir([inputDir filesep 'cAndH']);
for i=1:length(fileList)
    file = fileList(i).name;
    if ~strcmp(file(1),'.') && ~isempty(regexp(file,'FPKM.csv'))
        disp(file);
        FI = fopen([inputDir filesep 'cAndH' filesep file]);
        firstLine = fgetl(FI); words = strsplit(firstLine,','); observationIDs = words(2:end);
        dataFields = textscan(FI,['%s' repmat('%f',1,length(observationIDs))],'Delimiter',',','HeaderLines',1);
        fclose(FI);
        
        geneIDs = dataFields{1};
        expressionData = [ dataFields{2:end} ];
        
        [~, intersectIdxs, intersectIdxsAll] = intersect(observationIDs, cAndHObservationIDs);
        QLengths = cAndHQLengths(intersectIdxsAll);
        months = cAndHMonths(intersectIdxsAll);
        tissues = cAndHTissues(intersectIdxsAll);
        seqTypes = cAndHSeqTypes(intersectIdxsAll);
        expressionData = expressionData(:,intersectIdxs);
        cAndHExpressionData(:,intersectIdxsAll) = expressionData;
    end
end

cAndHGeneIDs = geneIDs;
save([outputDir1 filesep 'cAndHHDData.mat'],'cAndHObservationIDs', ...
    'cAndHTissues','cAndHQLengths','cAndHSexes','cAndHMonths','cAndHExpressionData','cAndHGeneIDs','cAndHSeqTypes');