config;

ensMito = {'ENSMUSG00000064341', ...
    'ENSMUSG00000064345', ...
    'ENSMUSG00000064360', ...
    'ENSMUSG00000064363', ...
    'ENSMUSG00000065947', ...
    'ENSMUSG00000064367', ...
    'ENSMUSG00000064368', ...
    'ENSMUSG00000064357', ...
    'ENSMUSG00000064356', ...
    'ENSMUSG00000064370', ...
    'ENSMUSG00000064351', ...
    'ENSMUSG00000064354', ...
    'ENSMUSG00000064358'};
commonNames = {'nd1','nd2','nd3','nd4','nd4l','nd5','nd6','atp6','atp8','cytb','co1','co2','co3'};
ensColors = linspace(1,10,length(ensMito));

count=0;
for i=1:length(uniqueMonths)
    for j=1:length(uniqueQLengths)
        count=count+1;
        xlabels{count} = ['Month' num2str(uniqueMonths(i)) 'QLength' num2str(uniqueQLengths(j))];
    end
end

for l=1:length(uniqueTissues)
    for k=1:length(ensMito)
        figure('Visible','off');
        count=0;
        vals = [];
        xvals = [];
        for i=1:length(uniqueMonths)
            for j=1:length(uniqueQLengths)
                count=count+1;
                disp([num2str(i) ' ' num2str(j)])
                [~, intersectIdxs1, intersectIdxs2] = intersect(allGeneIDs, ensMito{k});
                mask = allQLengths == uniqueQLengths(j) & allMonths == uniqueMonths(i) & strcmp(allSeqTypes,'mRNA');
                tempVals = allExpressionData(intersectIdxs1,mask);
                vals = [vals(:); tempVals(:)];
                xvals = [xvals(:); count*ones(length(tempVals(:)),1)];
            end
        end
        scatter(xvals(:),vals(:),[],'b');
        ylabel('RNASeq Counts');
        ylim([-.7*max(vals(:)) 1.3*max(vals(:))]);
        title(commonNames{k});
        for i=1:length(xlabels)
            text(i,-.7*max(vals(:)),xlabels{i},'Rotation',90);
        end
        saveas(gcf,[outputDir filesep uniqueTissues{l} filesep commonNames{k} '.png']);
        close(gcf);
    end
end