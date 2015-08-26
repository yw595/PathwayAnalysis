[~, ~, rawMasterDecoder] = xlsread( ...
['input' filesep 'Other Timepoints' filesep 'Master_Decoder_Ring.xlsx']);

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

if 1%~haveProcessedData
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
end

allGeneIDs = geneIDs;
save('AllHDData.mat','allObservationIDs','allMouseIDs','allTissues','allQLengths','allSexes','allMonths','allExpressionData','allGeneIDs','allSeqTypes');

uniqueTissues = {'cortex','liver','striatum'};
uniqueSexes = {'FEMALE','MALE'};
pValsCorr = zeros( 3,size(allExpressionData,1),length(uniqueTissues) );
pValsT = zeros( size(allExpressionData,1),length(uniqueTissues));
rhos = zeros( size(allExpressionData,1),length(uniqueTissues));

for i=1:size(allExpressionData,1)
    if mod(i,1000)==0
        disp(i);
    end
    count = 0;
    for k=1:length(uniqueTissues)
        count = count+1;
        mask = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA');
        dataRow = allExpressionData(i,:);
        mask1 = mask & allQLengths==20 & allMonths==2;
        mask2 = mask & allQLengths==175 & allMonths==10;
        [h pValsT(i,count)] = ttest2(dataRow(mask1)', dataRow(mask2)');
        
        mask1 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==2;
        [rho pValsCorr(1,i,count)] = corr(allQLengths(mask1), dataRow(mask1)', 'type', 'Spearman');
        mask2 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==6;
        [rho pValsCorr(2,i,count)] = corr(allQLengths(mask2), dataRow(mask2)', 'type', 'Spearman');
        mask3 = strcmp(allTissues, uniqueTissues{k}) & strcmp(allSeqTypes, 'mRNA') & allMonths==10;
        [rho pValsCorr(3,i,count)] = corr(allQLengths(mask3), dataRow(mask3)', 'type', 'Spearman');
    end
end

save('AllHDDataProcessed.mat','uniqueTissues','uniqueSexes','pValsCorr','pValsT','rhos');