function gapTable = sortGaps(rootGaps, model, param)
    indexes = [];
    ind = [];
    met_names = {};
    matches = {};
    if contains(param,'c0')
        matches = regexp(rootGaps, '\w*_c0\>', 'match' );
    elseif contains(param,'p')
        matches = regexp(rootGaps, '\w*_p\>', 'match' );
    elseif contains(param,'e0')
        matches = regexp(rootGaps, '\w*_e0\>', 'match' );
    elseif contains(param,'all')
        matches = regexp(rootGaps, '\<cpd\w*', 'match' );
    end
    matches = matches(~cellfun('isempty',matches));
    matches = cellstr(string(matches));
    for i=1:length(matches)
        k = matches(i);
        ind = find(contains(model.mets, k));
        indexes(i,1) = ind;
        met_names(i,1) = model.metNames(ind);
    end
gapTable = table(matches, met_names, indexes)
end