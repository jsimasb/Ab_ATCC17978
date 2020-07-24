% compareFluxes
% From one input model in two different conditions, find reactions that
% became blocked from condition one to condition two

function [FluxTable_AB, FluxTable_BA, nonZeroA, nonZeroB] = compareFluxes(modelA, modelB)
    if length(modelA.rxns) ~= length(modelB.rxns)
        display('Models must have the same set of reactions')
    end
    fluxA_AB = [];
    fluxB_AB = [];
    rxnID_AB = {};
    rxnEq_AB = {};
    
    fluxA_BA = [];
    fluxB_BA = [];
    rxnID_BA = {};
    rxnEq_BA = {};
    
    solutionA = solveLP(modelA);
    solutionB = solveLP(modelB);
    
    nonZeroA_ids = {};
    nonZeroA_eqs = {};
    nonZeroA_fluxes = [];
    
    nonZeroB_ids = {};
    nonZeroB_eqs = {};
    nonZeroB_fluxes = [];
    
    equations = constructEquations(modelA, modelA.rxns, false, false, false, true, false);
    % builds table AB
    for i=1:length(equations)
        if solutionA.x(i) ~= 0
            nonZeroA_ids(end+1,1) = modelA.rxns(i);
            nonZeroA_eqs(end+1,1) = equations(i);
            nonZeroA_fluxes(end+1,1) = solutionA.x(i);
            
            if solutionB.x(i) == 0
                fluxA_AB(end+1,1) = solutionA.x(i);
                fluxB_AB(end+1,1) = solutionB.x(i);
                rxnID_AB(end+1,1) = modelA.rxns(i);
                rxnEq_AB(end+1,1) = equations(i);
            end
        end
    end    
    FluxTable_AB = table(rxnID_AB, rxnEq_AB, fluxA_AB, fluxB_AB) ;
    nonZeroA = table(nonZeroA_ids, nonZeroA_eqs, nonZeroA_fluxes) ;
    
        % builds table BA
    for i=1:length(equations)
        if solutionB.x(i) ~= 0 
            nonZeroB_ids(end+1,1) = modelB.rxns(i);
            nonZeroB_eqs(end+1,1) = equations(i);
            nonZeroB_fluxes(end+1,1) = solutionB.x(i);
            
            if solutionA.x(i) == 0
                fluxA_BA(end+1,1) = solutionA.x(i);
                fluxB_BA(end+1,1) = solutionB.x(i);
                rxnID_BA(end+1,1) = modelA.rxns(i);
                rxnEq_BA(end+1,1) = equations(i);
            end
        end
    end
    FluxTable_BA = table(rxnID_BA, rxnEq_BA, fluxB_BA, fluxA_BA) ;
    nonZeroB = table(nonZeroB_ids, nonZeroB_eqs, nonZeroB_fluxes) ;
end