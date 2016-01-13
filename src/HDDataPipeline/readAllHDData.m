config;

[~, ~, rawMasterDecoder] = xlsread( ...
    [inputDir filesep 'Other Timepoints' filesep 'Master_Decoder_Ring.xlsx']);

% need to check format with 6 month data stored in HDData.mat
% in tissues, liver is lowercase here, QLengths do not have Q prefixed,
% sexes is fully written out
allObservationIDs = rawMasterDecoder(2:end,2);
allMouseIDs = rawMasterDecoder(2:end,1);
allTissues = rawMasterDecoder(2:end,3);
allQLengths = cell2mat(rawMasterDecoder(2:end,7));
allSexes = rawMasterDecoder(2:end,8);
allMonths = cell2mat(rawMasterDecoder(2:end,10));
allSeqTypes = rawMasterDecoder(2:end,4);

allExpressionData = NaN(23351,length(allObservationIDs));

fileList = dir(['input' filesep 'Other Timepoints']);
for i=1:length(fileList)
    file = fileList(i).name;
    % filter out ._ files, one near empty corresponding to each real
    % datafile, also use only FPKM, not rawcounts files
    if ~strcmp(file(1),'.') && ~isempty(regexp(file,'FPKM.xls'))
        %[junk1, ~, rawData] = xlsread( ...
        %['input' filesep 'Other Timepoints' filesep 'cortex_2m_FPKM.txt']);
        FI = fopen(['input' filesep 'Other Timepoints' filesep file]);
        firstLine = fgetl(FI); words = strsplit(firstLine,'\t'); observationIDs = words(4:end);
        fclose(FI);
        FI = fopen(['input' filesep 'Other Timepoints' filesep file]);
        dataFields = textscan(FI, ...
            ['%s%s%s' repmat('%f',1,length(observationIDs))],'Delimiter','\t','HeaderLines',1);
        fclose(FI);
        
        % use all arrays from master decoder above to get correct metadata
        % for this file
        geneIDs = dataFields{2};
        expressionDataHD = cell2mat(dataFields(4:(3+length(observationIDs))));
        
        [~, intersectIdxs, intersectIdxsAll] = intersect(observationIDs, allObservationIDs);
        observationIDs = observationIDs(intersectIdxs);
        mouseIDs = allMouseIDs(intersectIdxsAll);
        tissues = allTissues(intersectIdxsAll);
        QLengths = allQLengths(intersectIdxsAll);
        sexes = allSexes(intersectIdxsAll);
        months = allMonths(intersectIdxsAll);
        expressionDataHD = expressionDataHD(:,intersectIdxs);
        allExpressionData(:,intersectIdxsAll) = expressionDataHD;
    end
end

allGeneIDs = geneIDs;
save([inputDir filesep 'AllHDData.mat'],'allObservationIDs','allMouseIDs', ...
    'allTissues','allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');