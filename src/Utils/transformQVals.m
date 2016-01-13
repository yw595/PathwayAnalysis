function newQVals = transformQVals(oldQVals)
newQVals = zeros(size(oldQVals));
if analyzeAll
    map = [20 1;80 0;92 0;111 0;140 0;175 2];
else
    map = [20 1;80 2;92 3;111 4;140 5;175 6];
end
for index2=1:size(map,1)
    newQVals(oldQVals==map(index2,1))=map(index2,2);
end
end