configPicard;
v2struct();
disp('running readPicardData');
outputDir1 = [outputDir filesep 'readPicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

if isunix
    FI = fopen(['input' filesep 'Picard' filesep 'GSE56158_RPKM_cybrids_normalized_readgroups.txt']);
    firstLine = fgetl(FI); words = strsplit(firstLine,'\t'); PicardGeneIDs = words(2:end)';
    fclose(FI);
    FI = fopen(['input' filesep 'Picard' filesep 'GSE56158_RPKM_cybrids_normalized_readgroups.txt']);
    dataFields = textscan(FI, ...
        ['%s' repmat('%f',1,length(PicardGeneIDs))],'Delimiter','\t','HeaderLines',1);
    fclose(FI);
    PicardHeteroplasmies = dataFields{1};
    PicardHeteroplasmies = cellfun(@(x) readPicardDataHelper(x), PicardHeteroplasmies,'UniformOutput',0);
    PicardExpressionData = cell2mat(dataFields(2:end))';
end

save([outputDir1 filesep 'PicardData.mat'],'PicardHeteroplasmies','PicardExpressionData','PicardGeneIDs');