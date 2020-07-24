% validateModel

function [Out_struct, Out_table] = validateModel(model, cutoff)

% stores all files with txt extension in the current directory into Files struct
Files = dir("*.txt");

% initializes cell arrays CARBON_SOURCE and PHENOTYPE and array FLUX_BIOMASS
CARBON_SOURCE = {};
FLUX_BIOMASS = [];
PHENOTYPE = {};

% iterates through number of files with txt extension in the given directory
for k=1:length(Files)
  % stores current filename in a variable
  B = Files(k,1).name;

  % remove txt extension; creates condition name string
  condition = erase(B,'.txt');

  % stores growth condition being tested in a cell array of carbon sources
  CARBON_SOURCE(k,1) = cellstr(condition);

  % loads the current tab-delimited text file into a table representing the specific growth condition
  T = readtable(Files(k).name);

  % converts the second column in the table into a cell array
  C = table2cell(T(:,2));

  % adds _e0 to all elements of the array to create strings for compound IDs
  C= append(C, '_e0');

  % creates IDs for the reactions that will be added
  rxnIDs= append('EX_', C);
  
  % extracts compound names to a cell array of strings
  CPD_names = table2cell(T(:,1));
   
  % turns min flux column of table into array
  LB = table2array(T(:,4));
  
  % turns max flux column of table into array
  UB = table2array(T(:,5));

  % turns concentration column of table into array
  CONCENTRATION = table2array(T(:,3));
  
  add_rxnIDs = {};
  cpd_names = {};
  lb = [];
  ub = [];
  concentration = [];
  
  for i=1:length(rxnIDs)
    % filters compounds that are NOT already in current_model
    if contains(model.rxns, rxnIDs(i,1)) == 0
        add_rxnIDs(end+1) = rxnIDs(i,1);
        cpd_names(end+1) = C(i,1);
        lb(end+1) = LB(i,1);
        ub(end+1) = UB(i,1);
        concentration(end+1) = CONCENTRATION(i,1);
    end
  end
  
  % sets ids of reactions to be included
  rxnsToAdd.rxns = add_rxnIDs;

  % creates strings to name reactions
  %rxnNames= append('provides ', CPD_names, ' to the extracellular environment')

  % sets names of reactions
  rxnsToAdd.rxnNames = add_rxnIDs;

  % sets lower bounds of reactions to be included
  rxnsToAdd.lb = lb;

  % sets upper bounds of reactions to be included
  rxnsToAdd.ub = ub;

  % sets stoichiometric coefficients of reactions to be included
  %rxnsToAdd.stoichCoeffs = concentration   % NOT in use

  % creates equations that represent the sink reactions that will be added
  equations = strcat(concentration +" " + string(cpd_names) + " <=> ");

  % adds equations to input struct as a field
  rxnsToAdd.equations = cellstr(equations);
  
    % add the sink reactions to the model
    % model = addRxns(model, rxnsToAdd, 1, 'e0', false, false);
  new_model = addRxns(model, rxnsToAdd, 1, 'e0', true, false);

  % solves the linear programming problem of the model and stores flux values in the solution struct
  solution = solveLP(new_model);

  % stores the flux value through the objective function in a variable
  flux = solution.f;

  % stores flux through biomass reaction for each growth condition tested in an array
  FLUX_BIOMASS = [FLUX_BIOMASS; flux];

  % tests if the output represents Growth (G) or No Growth (NG) phenotypes and adds result for each condition in a cell array

  if abs(flux) >= cutoff
    PHENOTYPE{k,1} = 'G';
  else
    PHENOTYPE{k,1} = 'NG';
  end

% stores the main variables and structures of each condition of the Validation test in a struct
  Out_struct(k,1).conditionTested = condition;
  Out_struct(k,1).modelOutput = new_model;
  Out_struct(k,1).rxnInput = rxnsToAdd;
  Out_struct(k,1).solutionFound = solution;

end

% assembles simulation results in a table
Out_table = table(CARBON_SOURCE, FLUX_BIOMASS, PHENOTYPE);

% writes table to output text file
writetable(Out_table, 'Validation_result.csv')

end
