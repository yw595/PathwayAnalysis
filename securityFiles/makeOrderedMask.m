function orderedMaskArray = makeOrderedMask(mouseIDToInfo, mouseIDList, varargin)
    p = inputParser;
    p.FunctionName = 'orderedMaskArray';
    p.StructExpand = false;
    p.CaseSensitive = true;
    if isfield(p, 'PartialMatching')
        p.PartialMatching = false;
    end
    
    p.addParamValue('sex','');
    p.addParamValue('tissue','');
    p.addParamValue('QLength','');
    p.parse(varargin{:});
    results = p.Results;
    eval(structvars(numel(fieldnames(results)), results));
    
    orderedMaskArray = zeros(length(mouseIDList));
    for i=1:length(mouseIDList)
        mouseID = mouseIDList{i};
        if(strcmp(mouseIDToInfo(mouseID).tissue,tissue) && ...
            strcmp(mouseIDToInfo(mouseID).sex,sex) && ...
            mouseIDToInfo(mouseID).QLength==QLength)
            orderedMaskArray(i) = 1;
        end
    end
end