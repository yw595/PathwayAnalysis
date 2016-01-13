%configPicard;
%v2struct();
disp('running cytoscapePicardData');
outputDir1 = [outputDir filesep 'cytoscapePicardData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

FI = fopen([outputDir1 filesep 'recon2network.txt'],'w');
for i=1:length(origRecon2.rxns)
    metIdxs = find(origRecon2.S(:,i));
    for j=1:length(metIdxs)
        fprintf(FI,'%s\t%d\t%s\t%s',origRecon2.rxns{i},full(origRecon2.S(metIdxs(j),i)),origRecon2.mets{metIdxs(j)},origRecon2.subSystems{i});
        if i ~= length(origRecon2.rxns) || j~=length(metIdxs)
            fprintf(FI,'\n');
        end
    end
end
fclose(FI);