mouseEnsemblToHuman = table2cell(readtable('All Mouse Ensembl To Human.csv','Format','%s%s'));
humanEnsemblToEntrez = importdata('All Human Ensembl To Entrez.txt');
if ~exist('mouseEnsemblToEntrez','var')
mouseEnsemblToEntrez = containers.Map;
for i=1:size(mouseEnsemblToHuman,1)
    i
    lineNum = find(strcmp( humanEnsemblToEntrez.textdata(:,1),mouseEnsemblToHuman{i,2} ));
    if(length(lineNum)==1 && lineNum > 0)
        if(~isnan(humanEnsemblToEntrez.data(lineNum-1)))
            mouseEnsemblToEntrez(mouseEnsemblToHuman{i,1}) = humanEnsemblToEntrez.data(lineNum-1);
        end
    end
end
end

haveReadXLS = 1;
if exist('rawMetadata','var') && exist('rawData','var')
    haveReadXLS = 0;
end

if(haveReadXLS)
%get maps of mouseID to LineNum, QLength, Sex and Tissue
[junk1, junk2, rawMetadata] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA Metadata');
mouseIDToLineNum = containers.Map;
mouseIDToQLength = containers.Map;
mouseIDToSex = containers.Map;
mouseIDToTissue = containers.Map;
%TODO: see we can run cell and array index in one line
for i=2:size(rawMetadata,1)
    mouseIDToLineNum(rawMetadata{i,1}) = i;
    mouseIDToQLength(rawMetadata{i,1}) = str2num(rawMetadata{i,4}(2:end));
    mouseIDToSex(rawMetadata{i,1}) = rawMetadata{i,6};
    mouseIDToTissue(rawMetadata{i,1}) = rawMetadata{i,3};
end

[junk1, ~, rawData] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA FPKM data');
end

ensemblGeneIDArray = rawData(2:end,1);
hasEntrezID = isKey(mouseEnsemblToEntrez, ensemblGeneIDArray);
ensemblGeneIDArray = ensemblGeneIDArray(hasEntrezID);
entrezGeneIDArray = cellfun(@(x) mouseEnsemblToEntrez(x), ensemblGeneIDArray);

rawDataArray = cell2mat(rawData(2:end,2:end));
mouseIDArray = rawMetadata(2:end,1)';
tissueArray = rawMetadata(2:end,3)';
QLengthArray = rawMetadata(2:end,4)';
sexArray = rawMetadata(2:end,6)';

meanArray = zeros(size(rawDataArray,1),6); stdArray = zeros(size(rawDataArray,1),6);
tissues = {'cortex','Liver','striatum'};
sexes = {'F','M'};
QLengths = {'Q20','Q80','Q92','Q111','Q140','Q175'};
if ~exist('HDFalconData','dir')
    system('mkdir HDFalconData');
end
for i=1:length(tissues)
    for j=1:length(sexes)
        for k=1:length(QLengths)
            relevantData = rawDataArray(:, strcmp(tissueArray,tissues{i}) ...
            & strcmp(sexArray,sexes{j}) & strcmp(QLengthArray,QLengths{k}) );
            
            meanData = mean(relevantData,2);
            stdData = std(relevantData,0,2);
            meanData = meanData(hasEntrezID,:);
            stdData = stdData(hasEntrezID,:);
            
            outputFile = fopen(['HDFalconData' filesep 'HDFalconData' tissues{i} sexes{j} QLengths{k} '.txt'],'w');
            fprintf(outputFile,'entrez Gene ID\tmean\tstd\n');
            for l=1:length(entrezGeneIDArray)
                fprintf(outputFile,'%d\t%f\t%f',entrezGeneIDArray(l), meanData(l), stdData(l));
                if l~=length(entrezGeneIDArray)
                    fprintf(outputFile,'\n');
                end
            end
            fclose(outputFile);
        end
    end
end