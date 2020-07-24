% find all sink reactions in the model and iteratively remove one of them
% to identify which ones are essential for growth. 
matches = {};
% retrieves all sink reactions in the model
matches=regexp(current_modelo.rxns, '\<EX_\w*', 'match' );
% turns output in a cell array
matches = matches(~cellfun('isempty',matches));
matches = string(matches);
matches = cellstr(matches);
for i=1:length(matches)
    % Removes the sink reaction, without altering the input model
    new_model = removeRxns(current_modelo, matches(i), 'metFlag', false);
    % Runs FBA on new_model
    solution = solveLP(new_model);
    sinkAnalysis.flux(i,1) = solution.f;
    sinkAnalysis.sinkID(i,1) = matches(i);
    if solution.f == 0
        sinkAnalysis.class(i,1) = {'Essential'};
    else
        sinkAnalysis.class(i,1) = {'Non-essential'};
    end
end
    