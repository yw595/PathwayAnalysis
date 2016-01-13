function mediumExcIdxs = loadMediumExcIdxs2(model)
mediumExc = {
    'EX_ca2(e)','EX_so4(e)','EX_k(e)','EX_cl(e)','EX_na1(e)', ...
    'EX_hco3(e)','EX_pi(e)','EX_glc(e)','EX_gthrd(e)','EX_o2(e)', ...
    'EX_co2(e)','EX_h2o(e)','EX_h(e)'};
mediumExcIdxs = [];
for i = 1:length(mediumExc)
    mediumExcIdx = find(ismember(model.rxns, mediumExc{i}));
    if ~isempty(mediumExcIdx)
       mediumExcIdxs(end + 1) = mediumExcIdx;
    end
end
