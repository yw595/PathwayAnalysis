load('HDDataProcessed.mat');

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