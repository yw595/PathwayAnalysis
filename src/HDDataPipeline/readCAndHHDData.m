config;

[~, ~, rawCAndHMatrix] = xlsread([inputDir filesep 'GSE73508_series_matrix.txt']);
cAndHObservationIDs = rawCAndHMatrix(50,2:end);
rawQLengths = rawCAndHMatrix(30,2:end);
for i=1:length(rawQLengths)
    if isempty(regexp(rawQLengths{i},'WT'))
        x = rawQLengths{i};
        cAndHQLengths(i) = str2num(x(regexp(x,'(')+2:regexp(x,')')-1));
    else
        cAndHQLengths(i) = 0;
    end
end
cAndHMonths = cellfun(@(x) str2num(x(regexp(x,':')+2:regexp(x,'mon')-2)), rawCAndHMatrix(41,2:end));
cAndHTissues = cellfun(@(x) x(regexp(x,'-')+2:end), rawCAndHMatrix(43,2:end),'UniformOutput',0);
for i=2:size(rawCAndHMatrix,2)
    if ~isempty(regexp(rawCAndHMatrix{30,i},'mRNA'));
        cAndHSeqTypes{i-1} = 'mRNA';
    else
        cAndHSeqTypes{i-1} = 'miRNA';
    end
end
cAndHSexes = cellfun(@(x) upper(x(regexp(x,':')+2:end)), rawCAndHMatrix(42,2:end),'UniformOutput',0);

fileList = dir(['input' filesep 'Other Timepoints']);
for i=1:length(fileList)
    file = fileList(i).name;
    if ~strcmp(file(1),'.') && ~isempty(regexp(file,'FPKM.xls'))
        FI = fopen([inputDir filesep 'cAndH' filesep file]);
        firstLine = fgetl(FI); words = strsplit(firstLine,'\t'); observationIDs = words(2:end);
        dataFields = textscan(fi,['%s' repmat('%f',1,length(observationIDs))],'Delimiter',',','HeaderLines',1);
        fclose(fi);
        
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
save([inputDir filesep 'cAndHHDData.mat'],'cAndHObservationIDs', ...
    'cAndHTissues','cAndHQLengths','cAndHSexes','cAndHMonths','cAndHExpressionData','cAndHGeneIDs','cAndHSeqTypes');