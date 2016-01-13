config;
v2struct();
disp('running inputFALCONAllHHDData');
outputDir1 = [topDir filesep 'inputFALCONAllHHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

for i=1:length(uniqueTissues)
    for j=1:length(uniqueSexes)
        for k=1:length(uniqueQLengths)
            disp([uniqueTissues{i} ' ' uniqueSexes{j} ' ' uniqueQLengths{k}])
            relevantData = allGeneExpression(:, strcmp(allTissues,uniqueTissues{i}) ...
            & strcmp(allSexes,uniqueSexes{j}) & strcmp(allQLengths,uniqueQLengths{k}) );
            
            meanData = mean(relevantData,2);
            stdData = std(relevantData,0,2);
            
            outputFile = fopen(['HDFalconData' uniqueTissues{i} uniqueSexes{j} uniqueQLengths{k} '.txt'],'w');
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