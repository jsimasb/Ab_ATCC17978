
function [filteredByScore, alreadyInModel] = trimCandidates(candidate_reactions, model, gapScoreCutOff, metScoreCutOff)
    filteredByScore = struct([]);
    alreadyInModel = struct([]);

    alreadyInModel = isPresent(candidate_reactions, model);
    filteredByScore = filterScores(candidate_reactions,gapScoreCutOff, metScoreCutOff, model);

end

function filteredByScore = filterScores(candidate_reactions, gapScoreCutOff, metScoreCutOff, model)
    for i=1:length(candidate_reactions)
        filteredByScore(i,1).SEEDindex = [];
        filteredByScore(i,1).ids = {};
        filteredByScore(i,1).equations = {};
        filteredByScore(i,1).metScore = [];
        filteredByScore(i,1).gapScore = [];
        for h=1:length(candidate_reactions(i).ids)
            if candidate_reactions(i).metScore(h) >= metScoreCutOff && candidate_reactions(i).gapScore(h) >= gapScoreCutOff && isempty(find(contains(model.rxns, candidate_reactions(i).ids{h})))
                
                filteredByScore(i,1).SEEDindex(end+1,1) = candidate_reactions(i).SEEDindexes(h);
                filteredByScore(i,1).ids(end+1,1) = candidate_reactions(i).ids(h);
                filteredByScore(i,1).equations(end+1,1) = candidate_reactions(i).equations(h);
                filteredByScore(i,1).metScore(end+1,1) = candidate_reactions(i).metScore(h);
                filteredByScore(i,1).gapScore(end+1,1) = candidate_reactions(i).gapScore(h);
            end
        end
    end
end

function alreadyInModel = isPresent(candidate_reactions, model)
    for i=1:length(candidate_reactions)
      alreadyInModel(i,1).ids = {};
      alreadyInModel(i,1).equations = {};
      alreadyInModel(i,1).flux = [];
      alreadyInModel(i,1).gapScore = [];
      solution = solveLP(model);
      for h=1:length(candidate_reactions(i).ids)
        if ~isempty(find(contains(model.rxns, candidate_reactions(i).ids{h})))
          alreadyInModel(i,1).ids(end+1,1) = candidate_reactions(i).ids(h);
          alreadyInModel(i,1).equations(end+1,1) = candidate_reactions(i).equations(h);
          alreadyInModel(i,1).gapScore(end+1,1) = candidate_reactions(i).gapScore(h);
          % get flux for that reacion
          idx = find(contains(model.rxns, candidate_reactions(i).ids(h)));
          alreadyInModel(i,1).flux(end+1,1) = solution.x(idx);
        end
      end
    end
end
