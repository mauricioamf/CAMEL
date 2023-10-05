%
% Minimizes Etot - sum Ei*MW
%
%% Cleaning the workspace and the command window
clear;clc

%% Configurations for MATLAB
% changeCobraSolver('gurobi', 'LP');
% set(0,'DefaultFigureWindowStyle','docked')

%% Loading the enzyme-constrained model and general data
% load('../../1.Minimal_adjustment/Models/ecModels-master/ecYeastGEM/model/ecYeastGEM_batch.mat');
load('optStrain_kcat.mat')
modelSTR = optStrain;
clear optStrain

% load('../../1.Minimal_adjustment/Code_Yeast/Yu2021_std_010.mat')
% modelEXP = model;
% clear model
% 
load('enzUsageFVA.mat')
% 
% sp_title = 'Scatterplot of predicted and experimental values (log) - Yu2021\_std\_010';
fluxes_filename = 'fluxes_optStrain_kcat.csv';

%% Constraints for generating VS2 and ES2 (ETH)
% set media conditions
% cd heme_production_ecYeastGEM/code/GECKO/geckomat/kcat_sensitivity_analysis
% c_source = 'D-glucose exchange (reversible)';
% modelSTR = changeMedia_batch(modelSTR,c_source);
% cd ../../../../..

% fix specific growth rate at the dilution rate, allowing 1% flexibility
CS_MW = 0.18015;
Yield = 0.122;
V_bio = Yield*CS_MW;

modelSTR = changeRxnBounds(modelSTR, 'r_2111', V_bio*(1+0.01), 'u');
modelSTR = changeRxnBounds(modelSTR, 'r_2111', V_bio*(1-0.01), 'l');

% constrain target reaction
targetIndex = find(strcmpi(modelSTR.rxnNames,'heme exchange'));
WT_prod = 0.062;
tol = 1E-15;

modelSTR.lb(targetIndex) = (1-tol)*WT_prod;
modelSTR.ub(targetIndex) = (1+tol)*WT_prod;

ptotSTR = 0.5;
f = 0.5;
sigma = 0.5;

etotSTR = ptotSTR*f*sigma;

%% Construct the LP problem to minimize ES1
[~,nRxns] = size(modelSTR.S);
enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];
[enzymeNum,~] = size(enzymeIds);

LPproblem = buildLPproblemFromModel(modelSTR);

enzymes = modelSTR.enzymes(:);

pIdx = [];
rxnIdx = {};

for i=1:numel(enzymes)
    pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(enzymes(i))],"")));
    rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
end

% LPproblem.A(pIdx,enzymeIds) = LPproblem.A(pIdx,enzymeIds).*(etotSTR);
% LPproblem.ub(enzymeIds) = modelSTR.MWs;

LPproblem.c = zeros(size(LPproblem.A,2),1);
LPproblem.c(8144) = etotSTR;
LPproblem.c(enzymeIds) = -1;
LPproblem.c(enzymeIds) = modelSTR.MWs.*LPproblem.c(enzymeIds);

LPproblem.osense = 1;

%% Verify if the problem is a valid LP problem
fprintf('\n');
statusOK = verifyCobraProblem(LPproblem);
if statusOK == -1
    disp('Invalid LP problem')
end

%% Solve the problem
MinimizedFlux = solveCobraLP(LPproblem);

%% Get the solution(s)
enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_'))); 
enzymeIds(end,:) = [];

MinimizedFlux.x = MinimizedFlux.full;
MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds);

fprintf('\n');
fprintf('LP SOLUTION\n');
printFluxes(modelSTR, MinimizedFlux.x, true);

% fprintf('\n');
% fprintf('LP SOLUTION FOR ES \n');
% protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
% protNames.names = modelSTR.metNames(protNames.ids);
% printFluxes(modelSTR, MinimizedFlux.x, false, '', '', '', protNames.names);

%% Calculate correlation with experimental values
pred = {};
pred(:,1) = modelSTR.rxns(enzymeIds);
pred(:,1) = replace(pred(:,1), 'draw_prot_', '');
pred(:,2) = num2cell(MinimizedFlux.x(enzymeIds));
pred = cell2table(pred);
pred.Properties.VariableNames = {'Protein' 'Predicted'};
pred.Protein = char(pred.Protein);

enzUsageFVA.Protein = char(enzUsageFVA.Protein);

merged = innerjoin(pred, enzUsageFVA);
merged(ismember(merged.Predicted, 0),:)=[];
merged(ismember(merged.pUsage, 0),:)=[];

mergedLog = merged;
mergedLog.Predicted = log10(abs(mergedLog.Predicted));
mergedLog.pUsage = log10(abs(mergedLog.pUsage));

mergedLog.Predicted(isinf(mergedLog.Predicted)|isnan(mergedLog.Predicted)) = 0;
mergedLog.pUsage(isinf(mergedLog.pUsage)|isnan(mergedLog.pUsage)) = 0;

cor = {};
cor(1,:) = {'Pearson', 'Spearman', 'Kendall'};

[cor{2,1}, cor{3,1}] = corr(mergedLog.pUsage, mergedLog.Predicted, 'Type', char(cor(1,1)));
[cor{2,2}, cor{3,2}] = corr(mergedLog.pUsage, mergedLog.Predicted, 'Type', char(cor(1,2)));
[cor{2,3}, cor{3,3}] = corr(mergedLog.pUsage, mergedLog.Predicted, 'Type', char(cor(1,3)));

formatSpec = "%s's correlation: %f \t p-value: %s \n";
fprintf('\n');
fprintf('Correlation between CAMEL and original predictions:')
fprintf('\n');
fprintf(formatSpec, cor{:});

%% Calculate the RMSE
rmse = mergedLog;

rmse.RMSE = rmse.pUsage - rmse.Predicted;
rmse.RMSE = (rmse.RMSE).^2;
result = median(rmse.RMSE);
result = sqrt(result);

fprintf('\n');
formatSpecRMSE = "Log10 root median squared error: %f";
fprintf(formatSpecRMSE, result);
fprintf('\n');

%% How many proteins inside Es?
fprintf('\n');
[numE,~] = size(nonzeros(MinimizedFlux.x(enzymeIds)));
formatSpecNUM = "Number of predicted enzymes: %f";
fprintf(formatSpecNUM, numE);
fprintf('\n');

%% Plot scatterplot of predicted values vs. baseline
% figure
% hold on
% scatter(mergedLog.pUsage,mergedLog.Predicted, 'filled')
% % title(sp_title)
% xlabel('Previous values')
% ylabel('CAMEL-predicted values')
% p = polyfit(mergedLog.pUsage, mergedLog.Predicted, 1);
% px = [min(mergedLog.pUsage) max(mergedLog.pUsage)];
% py = polyval(p, px);
% R = corrcoef(mergedLog.pUsage, mergedLog.Predicted);
% str = join(['r = ', char(num2str(R(2)))], "");
% text(-5.5, -11.5, str);
% plot(px, py, 'LineWidth', 2);
% hold off

%% Get ratio
% protExp = cell2table(cell(0,0));
% protExp.Protein = modelEXP.enzymes;
% protExp.EsExp = modelEXP.ub(enzymeIds);
% protExp.EsExp(isinf(protExp.EsExp)) = NaN;
% protExp = rmmissing(protExp);
% protExp.Protein = char(protExp.Protein);
% 
% protPred = cell2table(cell(0,0));
% protPred.Protein = merged.Protein;
% protPred.Es = merged.pUsage(:);
% protPred.Properties.VariableNames = {'Protein' 'Es'};
% protPred.Protein = char(protPred.Protein);
% 
% mergedRatio = innerjoin(protExp, protPred);
% mergedRatio(ismember(mergedRatio.EsExp, 0),:)=[];
% mergedRatio(ismember(mergedRatio.Es, 0),:)=[];
% mergedRatio.ratio = mergedRatio.Es / mergedRatio.EsExp;

%% Get the metabolic flux of each reaction catalysed by the minimized enzymes
pred(ismember(pred.Predicted, 0),:)=[];

for i=1:size(pred.Protein(:,1))
    enzMin(i,1) = convertCharsToStrings(pred.Protein(i,:));
end

enzMinIdx = find(contains(modelSTR.enzymes, enzMin));
enzGenesMin = modelSTR.enzGenes(enzMinIdx);
% enzGenesMin = convertCharsToStrings(enzGenesMin);

modelSTRb = ravenCobraWrapper(modelSTR);
modelSTRb.grRules = modelSTRb.rules;
modelSTRb.geneNames = modelSTRb.genes;

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
RxnsFlux.Genes = char(modelSTRb.geneNames(enzMinIdx));
RxnsFlux.Protein = ProtMin;
RxnsFlux.Protein = char(RxnsFlux.Protein);
RxnsFlux.Fluxes = RxnsFluxArray;

%% Merge fluxes table with ratios table and export
mergedPred = innerjoin(pred, RxnsFlux);
writetable(mergedPred, fluxes_filename, 'Delimiter','\t')
fprintf('\n');
fprintf('Export finished');
fprintf('\n');