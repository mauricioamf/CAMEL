%
% Minimizes Etot - sum Ei*MW
%
%% Cleaning the workspace and the command window
clear;clc

%% Loading the enzyme-constrained model and general data
load('./Models_Ecoli/eciML1515_batch_<kcat/kapp>.mat');
modelSTR = ecModel_batch;
clear ecModel_batch

load('./Models_Ecoli/<growth condition>.mat')
modelEXP = ecModelP;
clear ecModelP

fluxes_filename = 'fluxes_<growth condition>_<kcat/kapp>.csv';

%% Constraints for generating VS2 and ES2 (ETH)
modelSTR = changeRxnBounds(modelSTR, 'EX_glc__D_e_REV', <value>, 'u');    % glucose

modelSTR = changeRxnBounds(modelSTR, 'EX_ac_e_REV', <value>, 'u');      % acetate
modelSTR = changeRxnBounds(modelSTR, 'EX_gam_e_REV', <value>, 'u');      % glucosamine
modelSTR = changeRxnBounds(modelSTR, 'EX_glyc_e_REV', <value>, 'u');       % glycerol
modelSTR = changeRxnBounds(modelSTR, 'EX_man_e_REV', <value>, 'u');       % mannose
modelSTR = changeRxnBounds(modelSTR, 'EX_pyr_e_REV', <value>, 'u');      % pyruvate
modelSTR = changeRxnBounds(modelSTR, 'EX_xyl__D_e_REV', <value>, 'u');       % xylose

modelSTR = FlexibilizeConstraints(modelSTR, <flexibization factor>);

% fix specific growth rate at the dilution rate, allowing 5% flexibility
modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', <value>*(1+0.5), 'u');
modelSTR = changeRxnBounds(modelSTR, 'BIOMASS_Ec_iML1515_core_75p37M', <value>*(1-0.5), 'l');

ptotSTR = <value>;

%% Predict enzyme usages
MinimizedFlux = PredictEnzymeUsage(modelSTR, ptotSTR);

%% Get ratio
enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];

protExp = cell2table(cell(0,0));
protExp.Protein = modelEXP.enzymes;
protExp.EsExp = modelEXP.ub(enzymeIds);
protExp.EsExp(isinf(protExp.EsExp)) = NaN;
protExp = rmmissing(protExp);
protExp.Protein = char(protExp.Protein);

protPred = cell2table(cell(0,0));
protPred.Protein = modelSTR.rxns(enzymeIds);
protPred.Protein = replace(protPred.Protein(:), 'draw_prot_', '');
protPred.Es = MinimizedFlux.x(enzymeIds);
protPred.Properties.VariableNames = {'Protein' 'Es'};
protPred.Protein = char(protPred.Protein);

mergedRatio = innerjoin(protExp, protPred);
%mergedRatio(ismember(mergedRatio.EsExp, 0),:)=[];
%mergedRatio(ismember(mergedRatio.Es, 0),:)=[];
mergedRatio.ratio = mergedRatio.Es ./ mergedRatio.EsExp;

%% Get the metabolic flux of each reaction catalysed by the minimized enzymes
for i=1:size(mergedRatio.Protein(:,1))
    enzMin(i,1) = convertCharsToStrings(mergedRatio.Protein(i,:));
end

enzMinIdx = find(contains(modelSTR.enzymes, enzMin));
enzGenesMin = modelSTR.enzGenes(enzMinIdx);

modelSTRb = ravenCobraWrapper(modelSTR);
rxnsEnz = findRxnsFromGenes(modelSTRb, enzGenesMin);

fieldNamesRxns = fieldnames(rxnsEnz);

for i=1:numel(fieldNamesRxns)
    for j=1:size(rxnsEnz.(fieldNamesRxns{i}))
        RxnNames(i,j) = rxnsEnz.(fieldNamesRxns{i})(j);
    end
end

for k=1:size(RxnNames(:,1))
    RxnsArray.(fieldNamesRxns{k}) = RxnNames(k,:);
    keep = any(~cellfun('isempty',RxnsArray.(fieldNamesRxns{k})), 1);
    RxnsArray.(fieldNamesRxns{k}) = RxnsArray.(fieldNamesRxns{k})(:,keep);
    RxnsArray.(fieldNamesRxns{k})(:,end) = [];
end

for y=1:size(RxnNames(:,1))
    RxnsIDArray.(fieldNamesRxns{y}) = find(contains(modelSTR.rxns, RxnsArray.(fieldNamesRxns{y})(:,:)));
end

for z=1:size(RxnNames(:,1))
    RxnsFluxArray.(fieldNamesRxns{z}) = MinimizedFlux.x(RxnsIDArray.(fieldNamesRxns{z})(:,:));
    RxnsFluxArray.(fieldNamesRxns{z}) = nonzeros(RxnsFluxArray.(fieldNamesRxns{z}));
    RxnsFluxArray.(fieldNamesRxns{z}) = max(RxnsFluxArray.(fieldNamesRxns{z}));
end

RxnsFluxArray = struct2cell(RxnsFluxArray);
ProtMin = modelSTR.enzymes(enzMinIdx);

RxnsFlux = cell2table(cell(0,0));
RxnsFlux.Protein = ProtMin;
RxnsFlux.Fluxes = RxnsFluxArray;
RxnsFlux.Protein = char(RxnsFlux.Protein);

%% Merge fluxes table with ratios table and export
mergedRatio = innerjoin(mergedRatio, RxnsFlux);
writetable(mergedRatio, fluxes_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');
