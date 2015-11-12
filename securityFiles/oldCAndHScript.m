% fi = fopen([inputDir filesep 'GSE73468_mRNA_2__6__and_10_month_cerebellum_FPKM.csv']);
% headerCere = fgetl(fi); headerCere = headerCere(2:end);
% rawCere = textscan(fi,['%s' repmat('%f',1,166)],'Delimiter',',','HeaderLines',1);
% fclose(fi);
% fi = fopen([inputdir filesep 'input\GSE73503_mRNA_2__6__and_10_month_hippocampus_FPKM.csv']);
% headerHippo = fgetl(fi); headerHippo = headerHippo(2:end);
% rawHippo = textscan(fi,['%s' repmat('%f',1,168)],'Delimiter',',','HeaderLines',1);
% fclose(fi);
% 
% [~, ~, intersectIdxs2] = intersect(cAndHObservationIDs, headerHippo);
% hippoQLengths = cAndHQLengths(intersectIdxs2);
% hippoMonths = cAndHMonths(intersectIdxs2);
% hippoTissues = cAndHTissues(intersectIdxs2);
% hippomRNA = cAndHmRNA(intersectIdxs2);
% [~, ~, intersectIdxs2] = intersect(cAndHObservationIDs, headerCere);
% cereQLengths = cAndHQLengths(intersectIdxs2);
% cereMonths = cAndHMonths(intersectIdxs2);
% cereTissues = cAndHTissues(intersectIdxs2);
% ceremRNA = cAndHmRNA(intersexctIdxs2);
%
%dataCere = [ rawCere{2:end} ];
%dataHippo = [ rawHippo{2:end} ];