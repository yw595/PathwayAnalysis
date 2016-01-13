inputFile = fopen('hLGDB Default Original.csv');

entrezIDs = [];
line = fgetl(inputFile);
line = fgetl(inputFile); % skip header of inputFile
while line ~= -1
    tf = isletter(line);
    %skip weird entry at line 189, 
    %some character with ASCII 9 at beginning
    if(~strcmp(line(2:end),'21752829'))
        entrezIDs(end+1) = str2num(line(1:find(tf,1)-1));
    end
    line = fgetl(inputFile);
end
fclose(inputFile);

outputFile = fopen('hLGDB Default Human Entrez List.txt','w');
for i=1:length(entrezIDs)-1
    fprintf(outputFile,'%d\n',entrezIDs(i));
end
fprintf(outputFile,'%d',entrezIDs(i));
fclose(outputFile);