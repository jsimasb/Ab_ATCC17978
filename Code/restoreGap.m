% restoreGap
% From input set of gaps, add a sink reaction for each gap (iteratively)
% Run FBA and check if the addition of the sink reaction altered model
% functionality

function Restored_table = restoreGap(model, input_gaps)
    restored_gap = {};
    new_flux = [];
    for i=1:length(input_gaps)
        rxnID = append('EX_', input_gaps(i));
        rxnsToAdd.rxns = rxnID;
        rxnsToAdd.rxnNames = rxnID;

        rxnEquation = append(input_gaps(i), ' <=> ');
        rxnsToAdd.equations = rxnEquation;

        new_model = addRxns(model, rxnsToAdd, 1);

        solution = solveLP(new_model);

        restored_gap(i,1) = input_gaps(i);
        new_flux(i,1) = solution.f;
    end
    
    Restored_table = table(restored_gap, new_flux);
end
    