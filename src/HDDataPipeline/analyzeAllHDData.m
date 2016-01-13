config;
v2struct();
disp('running analyzeAllHDData');
outputDir1 = [topDir filesep 'analyzeAllHDData'];
if ~exist(outputDir1,'dir')
    mkdir(outputDir1)
end

testT = 0; testCorr = 1;

if(testT)
    pVals = pValsT;
    labels = labelsT;
end
if(testCorr)
    pVals = pValsCorr;
    labels = labelsCorr;
end

pVals(isnan(pVals)) = 1;
sizeLabels = cellfun(@(x) length(x), labels);
labelIdxs = ones(length(labels),1);
while labelIdxs ~= -1
    labelIdxs
    mask = zeros(size(pVals));
    command = ('mask(:');
    for i=1:length(labelIdxs)
        command = [command ',' num2str(labelIdxs(i))];
    end
    command = [command ') = 1;'];
    eval(command); mask = (mask~=0);
    
    [sortPValsCol idxs] = sort(pVals(mask));
    cutoff = FDRCutoff(sortPValsCol);
    sigEntrezIDs = allGeneIDs(idxs(1:cutoff));
    sortPValsCol = sortPValsCol(1:cutoff);
    filename = 'sigEntrezIDs'; 
    for i=1:length(labelIdxs)
        ithLabels = labels{i};
        if iscell(labels{i})
            label = ithLabels{labelIdxs(i)};
        else
            label = num2str(ithLabels(i));
        end
        filename = [filename label]; 
    end
    writeData({sigEntrezIDs, sortPValsCol},[outputDir1 filename '.txt'],'\t');
    labelIdxs = makePermutes(sizeLabels,labelIdxs);
end
% for i=1:size(pVals,2)
%     fid = fopen([outputDir1 filesep 'sigEntrezIDs' num2str(i) '.txt'],'w');
%     ithSigEntrezIDs = sigEntrezIDs{i};
%     for j=1:length(ithSigEntrezIDs)
%         fprintf(fid,'%s\n',ithSigEntrezIDs{j});
%     end
%     fclose(fid);
% end
% 
% for k=1:length(uniqueTissues)
%     ithSigEntrezIDs = sigEntrezIDs{k};
%     if testCorr
%         mask = strcmp(allTissues, uniqueTissues{k}) ...
%             & allMonths==allMonths(timePoint) & strcmp(allSeqTypes,'mRNA');
%     else
%         mask = strcmp(allTissues, uniqueTissues{k}) & (allQLengths==20 | allQLengths==175) ...
%             & (allMonths==2 | allMonths==10) & strcmp(allSeqTypes,'mRNA');
%     end
%     [~, intersectIdxs, ~] = intersect(ithLookSigEntrezIDs,allGeneIDs);
%     sigEntrezVals{k} = allExpressionData( intersectIdxs,mask );
%     sigEntrezQLengths{k} = allQLengths(mask);
% end
% 
% for i=1:length(lookSigEntrezQLengths)
%     lookSigEntrezQLengths{i} = transformQVals(lookSigEntrezQLengths{i});
% end
% 
% save([outputDir1 filesep 'analyzeAllHDData.mat'], 'sigEntrezQLengths', 'sigEntrezVals', 'sigEntrezIDs', 'sortPValsColArray');