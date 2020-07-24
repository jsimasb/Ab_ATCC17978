% findSEEDCandidateRxns

function candidate_reactions = findCandidateRxns(model, gaps, allGaps, database)
  gaps = cellstr(gaps);
  Index = [];
  % Loop through input gaps array
  for k=1:length(gaps)
    % get the gap metabolite id
    gapmet = gaps(k);
    candidate_reactions_count = 0;

    % Find indexes in the SEED DATABASE for the gap metabolite
    Index = find(contains(database.equations, gapmet));

    % checks if candidate reactions were found
    if ~isempty(Index)
      % Loop through all candidate reactions found
      for i=Index
        % store candidate reaction names in cell array
        candidate_reactions_names = database.rxnNames(i);
        % store candidate reaction ids in array
        candidate_reactions_ids = database.rxns(i);
        % stores candidate reaction equations in array
        candidate_reactions_equations =  database.equations(i);
      end
      % Store the root gap involved in the candidate reactions
      candidate_reactions(k,1).solves = gaps(k);
      % Store reaction indexes in output structure
      candidate_reactions(k,1).SEEDindexes = Index;
      % Store candidate reaction count in output struct
      candidate_reactions(k,1).count = length(candidate_reactions(k,1).SEEDindexes);
      % Store candidate reaction names in output struct
      candidate_reactions(k,1).names = candidate_reactions_names;
      % Store candidate reaction identifiers in output struct
      candidate_reactions(k,1).ids = candidate_reactions_ids;
      % Store candidate reaction equations in output struct
      candidate_reactions(k,1).equations = candidate_reactions_equations;

    % If no candidate reactions were found
    else
      candidate_reactions(k,1).solves = {};
      candidate_reactions(k,1).SEEDindexes = [];
      candidate_reactions(k,1).count = [];
      candidate_reactions(k,1).names = {};
      candidate_reactions(k,1).ids = {};
      candidate_reactions(k,1).equations = {};
    end
   end

   for j=1:length(candidate_reactions)
     % If at least one candidade reaction was found
     if ~isempty(candidate_reactions(j).ids)
         % obtains all compounds in the reaction equation
         compounds=regexp(candidate_reactions(j).equations, '\<cpd\w*', 'match' );
         if ~isempty(compounds)
             for m=1:length(compounds)
               % Reset counters
               metcount = 0;
               notGapCount = 0;
               % checks if compound is present in the input model - Criterion 1
               for n=1:length(compounds{m})
                   is_present = find(contains(model.mets, compounds{m}(n)));
                   % if compound already exists in model, add 1 to metcount
                   if ~isempty(is_present)
                     metcount = metcount +1;
                     candidate_reactions(j).metToAdd(m,n) = {''};
                     % checks if compound in the reaction equation is a gap metabolite - Criterion 2
                     is_gap = find(contains(allGaps, compounds{m}(n)));
                     % if metabolite is NOT a gap in the model, this means, if it is connected to the central network
                         if isempty(is_gap)
                             notGapCount = notGapCount + 1;
                         end
                    % if compound does not exist in the input model, add it to the metToAdd field
                     else
                     candidate_reactions(j).metToAdd(m,n) = compounds{m}(n);
                     end
               end
               % rate of metabolites already in input model per total metabolites in equation, excluding the root gap
               Score.metScore(m,1) = (metcount - 1)/(length(compounds{m}) - 1);
               % rate of metabolites connected to the central network per total metabolites in equation
               Score.gapScore(m,1) = (notGapCount)/(length(compounds{m}) -1 );
               
               candidate_reactions(j).metScore(m,1) = Score.metScore(m,1);
               candidate_reactions(j).gapScore(m,1) = Score.gapScore(m,1);
             end
         end
           else
              candidate_reactions(j).metScore(m,1) = [];
              candidate_reactions(j).gapScore(m,1) = [];
              candidate_reactions(j).metToAdd(m,1) = {};
           end
     end
end
