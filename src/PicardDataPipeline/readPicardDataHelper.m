function y = readPicardDataHelper(x)

y = strrep(x,' ','');
if strcmp(x(end),'%') 
    y=x(1:end-1);
end
end