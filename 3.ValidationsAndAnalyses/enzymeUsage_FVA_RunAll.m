%
% Minimizes Etot - sum Ei*MW
%
%% Cleaning the workspace and the command window
clear;clc

%% Configurations for MATLAB
% initCobraToolbox(false);
changeCobraSolver('gurobi', 'LP');
%run("~/Softwares/RAVEN/installation/checkInstallation.m");
warning("off");

clear;clc

%% Loop over growth conditions for ecYeast8 kcat model

fprintf('\n' + "Starting FVA for kcat paremeterized ecYeast8 batch model");

% Constraints for suboptimal conditions
load('./constraints_Scerevisiae2.mat')

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:end);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('../../1.Minimal_adjustment/Models/ecModels-master/ecYeastGEM/model/ecYeastGEM_batch.mat');               
    modelSTR = ecModel_batch;
    clear ecModel_batch
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.05);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(13);

    f = 0.5;
    sigma = 0.5;
    etotSTR = ptotSTR*f*sigma;

    % modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Perform Flux Variability Analysis on enzyme usage
    FVAtable = enzymeUsage_FVA(modelSTR);

    % Export results table
    FVA_filename = "FVA_yeast_kcat_" + current_condition + ".csv";
    writetable(FVAtable, FVA_filename, 'Delimiter','\t')
    fprintf('\n');
    fprintf('Export finished');
    fprintf('\n');

end

%% Loop over growth conditions for ecYeast8 kapp model

clear;

fprintf('\n' + "Starting FVA for kapp paremeterized ecYeast8 batch model");

% Constraints for suboptimal conditions
load('./constraints_Scerevisiae2.mat')

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:end);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('../../1.Minimal_adjustment/Code_Yeast/ecYeastGEM_batch_kapp.mat')                
    modelSTR = ecModel_batch;
    clear ecModel_batch
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'r_1672', flux_values(3), 'u');       % CO2
    modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', flux_values(4), 'u');   % O2
    modelSTR = changeRxnBounds(modelSTR, 'r_2033', flux_values(5), 'u');      % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'r_2056', flux_values(6), 'u');      % succinate
    modelSTR = changeRxnBounds(modelSTR, 'r_1808', flux_values(7), 'u');       % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'r_1634', flux_values(8), 'u');      % acetate
    modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', flux_values(9), 'u');   % ethanol

    modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', flux_values(10), 'u');       % Gln
    modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', flux_values(11), 'u');       % Phe
    modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', flux_values(12), 'u');       % Ile

    modelSTR = FlexibilizeConstraints(modelSTR, 0.05);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1+0.05), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'r_2111', flux_values(1)*(1-0.05), 'l');

    ptotSTR = flux_values(13);

    f = 0.5;
    sigma = 0.5;
    etotSTR = ptotSTR*f*sigma;

    % modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Perform Flux Variability Analysis on enzyme usage
    FVAtable = enzymeUsage_FVA(modelSTR);

    % Export results table
    FVA_filename = "FVA_yeast_kapp_" + current_condition + ".csv";
    writetable(FVAtable, FVA_filename, 'Delimiter','\t')
    fprintf('\n');
    fprintf('Export finished');
    fprintf('\n');

end

%% Loop over growth conditions for eciML1515 kcat model

clear;

fprintf('\n' + "Starting FVA for kcat paremeterized eciML1515 batch model");

% Constraints for suboptimal conditions
load('./constraints_Ecoli.mat')

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:end);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('../../1.Minimal_adjustment/Code_Ecoli/eciML1515_batch.mat')                
    modelSTR = ecModel_batch;
    clear ecModel_batch
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.5), 'l');
    
    ptotSTR = flux_values(9);

    f = 0.5;
    sigma = 0.5;
    etotSTR = ptotSTR*f*sigma;

    % modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Perform Flux Variability Analysis on enzyme usage
    FVAtable = enzymeUsage_FVA(modelSTR);

    % Export results table
    FVA_filename = "FVA_Ecoli_kcat_" + current_condition + ".csv";
    writetable(FVAtable, FVA_filename, 'Delimiter','\t')
    fprintf('\n');
    fprintf('Export finished');
    fprintf('\n');

end

%% Loop over growth conditions for eciML1515 kapp model

clear;

fprintf('\n' + "Starting FVA for kapp paremeterized eciML1515 batch model");

% Constraints for suboptimal conditions
load('./constraints_Ecoli.mat')

% Define the names of the growth conditions
growth_conditions = constraints.Properties.VariableNames(3:end);

for k = 1:numel(growth_conditions)
        
    % Get the current growth condition
    current_condition = growth_conditions{k};

    fprintf('\n' + "Running simulations for " + current_condition);

    % Get the flux values for the current growth condition
    flux_values = constraints{:, current_condition};

    % pcGEM without integrated proteomics data
    load('../../1.Minimal_adjustment/Code_Ecoli/eciML1515_batch_kapp.mat')                
    modelSTR = ecModel_batch;
    clear ecModel_batch
    
    % Constraints for generating VS2
    modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', flux_values(2), 'u');    % glucose
    modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', flux_values(3), 'u');       % acetate
    modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', flux_values(4), 'u');   % glucosamine
    modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', flux_values(5), 'u');      % glycerol
    modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', flux_values(6), 'u');      % mannose
    modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', flux_values(7), 'u');       % pyruvate
    modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', flux_values(8), 'u');      % xylose

    modelSTR = FlexibilizeConstraints(modelSTR, 0.9);

    % fix specific growth rate at the dilution rate, allowing 5% flexibility
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1+0.5), 'u');
    modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', flux_values(1)*(1-0.5), 'l');
    
    ptotSTR = flux_values(9);

    f = 0.5;
    sigma = 0.5;
    etotSTR = ptotSTR*f*sigma;

    % modelSTR = changeRxnBounds(modelSTR, 'prot_pool_exchange', etotSTR, 'b');
    
    % Perform Flux Variability Analysis on enzyme usage
    FVAtable = enzymeUsage_FVA(modelSTR);

    % Export results table
    FVA_filename = "FVA_Ecoli_kapp_" + current_condition + ".csv";
    writetable(FVAtable, FVA_filename, 'Delimiter','\t')
    fprintf('\n');
    fprintf('Export finished');
    fprintf('\n');

end
