config;

for i=1:length(uniqueTissues)
    for j=1:length(uniqueSexes)
        for k=1:length(uniqueQLengths)
            relevantData = allGeneExpression(:, strcmp(tissueArray,tissues{i}) ...
            & strcmp(sexArray,sexes{j}) & strcmp(QLengthArray,QLengths{k}) );
            
            meanData = mean(relevantData,2);
            stdData = std(relevantData,0,2);
            
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