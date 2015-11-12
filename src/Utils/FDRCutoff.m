function cutoff = FDRCutoff(sortPvals)
cutoff = 1;
for index1=1:length(sortPvals)
    if(sortPvals(index1)<=index1*.05/length(sortPvals))
        cutoff=cutoff+1;
    else
        break;
    end
end
end