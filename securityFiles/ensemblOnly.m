inputFile = fopen('Ensembl Default List.txt');

ensemblIDs = {};
line = fgetl(inputFile);
line = fgetl(inputFile); % skip header of inputFile
while line ~= -1
    words = strsplitYiping(line,'\t');
    possibleEnsemblID = words{2};
    %check if actual ENSG, might be LRG
    if(strcmp(possibleEnsemblID(1:3),'ENS'))
        ensemblIDs{end+1} = possibleEnsemblID;
    end
    line = fgetl(inputFile);
end
fclose(inputFile);

outputFile = fopen('Ensembl (Only) Default List.txt','w');
for i=1:length(ensemblIDs)-1
    fprintf(outputFile,'%s\n',ensemblIDs{i});
end
fprintf(outputFile,'%s',ensemblIDs{i});
fclose(outputFile);