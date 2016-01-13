config;
v2struct();
disp('running readThreeHDData');
outputDir1 = [outputDir filesep 'readThreeHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if isunix
    FI = fopen(['input' filesep 'Other Timepoints' filesep 'Master_Decoder_Ring.txt']);
    dataFields = textscan(FI, ...
        '%s%s%s%s%s%s%s%s%s%s%s','Delimiter','\t','HeaderLines',1);
    fclose(FI);
    rawMasterDecoder = [dataFields{:}];
    threeObservationIDs = rawMasterDecoder(1:end,2);
    threeMouseIDs = rawMasterDecoder(1:end,1);
    threeTissues = rawMasterDecoder(1:end,3);
    threeQLengths = cellfun(@(x) [x 'Q'], rawMasterDecoder(1:end,7),'UniformOutput',0);
    threeSexes = rawMasterDecoder(1:end,8);
    threeMonths = cellfun(@(x) [x 'Months'], rawMasterDecoder(1:end,10),'UniformOutput',0);
    threeSeqTypes = rawMasterDecoder(1:end,4);
else
    [~, ~, rawMasterDecoder] = xlsread( ...
        [inputDir filesep 'Other Timepoints' filesep 'Master_Decoder_Ring.txt']);
    threeObservationIDs = rawMasterDecoder(2:end,2);
    threeMouseIDs = rawMasterDecoder(2:end,1);
    threeTissues = rawMasterDecoder(2:end,3);
    threeQLengths = cellfun(@(x) [x 'Q'], rawMasterDecoder(2:end,7));
    threeSexes = rawMasterDecoder(2:end,8);
    threeMonths = cellfun(@(x) [x 'Months'], rawMasterDecoder(2:end,10));
    threeSeqTypes = rawMasterDecoder(2:end,4);
end

threeExpressionData = NaN(23351,length(threeObservationIDs));

fileList = dir([inputDir filesep 'Other Timepoints']);
for i=1:length(fileList)
    file = fileList(i).name;
    if ~strcmp(file(1),'.') && ~isempty(regexp(file,'FPKM.xls'))
        disp(file);
        FI = fopen(['input' filesep 'Other Timepoints' filesep file]);
        firstLine = fgetl(FI); words = strsplit(firstLine,'\t'); observationIDs = words(4:end);
        fclose(FI);
        FI = fopen(['input' filesep 'Other Timepoints' filesep file]);
        dataFields = textscan(FI, ...
            ['%s%s%s' repmat('%f',1,length(observationIDs))],'Delimiter','\t','HeaderLines',1);
        fclose(FI);
        
        % use three arrays from master decoder above to get correct metadata
        % for this file
        geneIDs = dataFields{2};
        expressionDataHD = cell2mat(dataFields(4:(3+length(observationIDs))));
        
        [~, intersectIdxs, intersectIdxsThree] = intersect(observationIDs, threeObservationIDs);
        observationIDs = observationIDs(intersectIdxs);
        mouseIDs = threeMouseIDs(intersectIdxsThree);
        tissues = threeTissues(intersectIdxsThree);
        QLengths = threeQLengths(intersectIdxsThree);
        sexes = threeSexes(intersectIdxsThree);
        months = threeMonths(intersectIdxsThree);
        expressionDataHD = expressionDataHD(:,intersectIdxs);
        threeExpressionData(:,intersectIdxsThree) = expressionDataHD;
    end
end

threeGeneIDs = geneIDs;
save([outputDir1 filesep 'threeHDData.mat'],'threeObservationIDs','threeMouseIDs', ...
    'threeTissues','threeQLengths','threeSexes','threeMonths','threeExpressionData','threeGeneIDs','threeSeqTypes');