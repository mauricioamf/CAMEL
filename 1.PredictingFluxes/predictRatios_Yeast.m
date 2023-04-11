%
% Minimizes Etot - sum Ei*MW
%
%% Cleaning the workspace and the command window
clear;clc

%% Loading the enzyme-constrained model and general data
load('./Models_Scerevisiae/ecYeastGEM_batch_<kcat/kapp>.mat')
modelSTR = ecModel_batch;
clear ecModel_batch

load('./Models_Scerevisiae/<growth_condition>.mat')
modelEXP = model;
clear model

fluxes_filename = 'fluxes_<growth_condition>_<kcat/kapp>.csv';

%% Constraints for generating VS2 and ES2 (ETH)
modelSTR = changeRxnBounds(modelSTR, 'r_1714_REV', <value>, 'u');    % glucose
modelSTR = changeRxnBounds(modelSTR, 'r_1672', <value>, 'u');       % CO2
modelSTR = changeRxnBounds(modelSTR, 'r_1992_REV', <value>, 'u');   % O2
modelSTR = changeRxnBounds(modelSTR, 'r_2033', <value>, 'u');      % pyruvate
modelSTR = changeRxnBounds(modelSTR, 'r_2056', <value>, 'u');      % succinate
modelSTR = changeRxnBounds(modelSTR, 'r_1808', <value>, 'u');       % glycerol
modelSTR = changeRxnBounds(modelSTR, 'r_1634', <value>, 'u');      % acetate
modelSTR = changeRxnBounds(modelSTR, 'r_1761_REV', <value>, 'u');   % ethanol

modelSTR = changeRxnBounds(modelSTR, 'r_1891_REV', <value>, 'u');       % Gln
modelSTR = changeRxnBounds(modelSTR, 'r_1903_REV', <value>, 'u');       % Phe
modelSTR = changeRxnBounds(modelSTR, 'r_1897_REV', <value>, 'u');       % Ile

modelSTR = FlexibilizeConstraints(modelSTR, <flexibization factor);

% fix specific growth rate at the dilution rate, allowing 1% flexibility
modelSTR = changeRxnBounds(modelSTR, 'r_2111', <value>*(1+0.01), 'u');
modelSTR = changeRxnBounds(modelSTR, 'r_2111', <value>*(1-0.01), 'l');

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