haveReadXLS = 1;
if exist('rawMetadata','var') && exist('rawData','var')
    haveReadXLS = 0;
end

if(haveReadXLS)
    [junk1, junk2, rawMetadata] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA Metadata');
    mouseIDs = rawMetadata(2:end,1);
    tissues = rawMetadata(2:end,3);
    QLengths = rawMetadata(2:end,4);
    sexes = rawMetadata(2:end,6);

    [junk1, ~, rawData] = xlsread( ...
    '6_month_allelic_series_data_mRNA_with_metadata2014_21.xlsx','mRNA FPKM data');
    geneIDs = rawData(2:end,1);
end
expressionData = cell2mat(rawData(2:end,2:end));

GCTFile = fopen('HD.gct','w');
fprintf(GCTFile,'#1.2\n');
fprintf(GCTFile,'1000\t130\n');
fprintf(GCTFile,'NAME\tDescription\t');
for i=1:length(mouseIDs)
    if i==length(mouseIDs)
        fprintf(GCTFile,'%s\n',mouseIDs{i});
    else
        fprintf(GCTFile,'%s\t',mouseIDs{i});
    end
end
for i=1:length(geneIDs)
    fprintf(GCTFile,'%s\t',geneIDs{i});
    fprintf(GCTFile,'na\t');
    for j=1:size(expressionData,2)
        if j==size(expressionData,2)
            fprintf(GCTFile,'%f\n',expressionData(i,j));
        else
            fprintf(GCTFile,'%f\t',expressionData(i,j));
        end
    end
end
fclose(GCTFile);