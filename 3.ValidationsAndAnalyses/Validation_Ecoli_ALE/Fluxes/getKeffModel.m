%% Cleaning the workspace and the command window
clear;clc
initCobraToolbox(false);
tic
%% Loading the enzyme-constrained model and general data
% load('../../1.Minimal_adjustment/Code_Ecoli/eciML1515_batch.mat');

%% Generating keff adjusted model
% current = pwd;
% cd ../../1.Minimal_adjustment/PRESTO/Program
% 
% ecModel_batch.name = 'Escherichia coli';
% enzMetPfx = 'prot_';
% 
% ecModel_batch_keff = create_pFBAmod(ecModel_batch, enzMetPfx);
% 
% cd (current)

%% Generates keff-adjusted model and predicts enzyme usage
load('./sampleIDs.mat');

resultsGEM = cell2table(cell(0,0));

for k = 1:numel(sampleIDs)

    clearvars -except sampleIDs k resultsGEM

    load('../../1.Minimal_adjustment/Code_Ecoli/eciML1515_batch.mat');

    currentSample = sampleIDs{k};
    
    fprintf('\n' + "Building keff-adjusted model for  " + currentSample + '\n');
    
    raw_batchmod = ecModel_batch;
    enzMetPfx = 'prot_';
    
    kappFilename = "data_" + currentSample + ".xlsx";
    
    kapp=readtable(kappFilename, 'Sheet', 'Sheet2');
    %adapt naming
    kapp.Properties.VariableNames([1 4 5])={'ReactionID', 'kappmax', 'abundance'};
    %get reaction mapping information
    map=readtable(kappFilename, 'Sheet', 'Sheet3');
    %filter unmatched
    map=map(~cellfun(@isempty, map.final_match),[1,4]);
    %duplicate isozyme entries. 
    
    i=1;
    while i<=size(map,1)
        if contains(map.final_match(i), ';')
            isozms=strsplit(map.final_match{i}, '; ')';
            tmp_map=map(repelem(i, length(isozms)), :);
            tmp_map.final_match=isozms;
            map=[map(1:(i-1), :); tmp_map; map((i+1):end,:)];
            i=i+length(isozms);
        else
            i=i+1;
        end
    end
    
    enzMetIdx = find(contains(raw_batchmod.mets, enzMetPfx) &~ismember(raw_batchmod.mets, {'prot_pool'}));
    
    davidi_mod=raw_batchmod;
    
    for i=1:size(map,1)
        %retrieve davidi kapp
        new_kapp=-1/(kapp.kappmax(ismember(kapp.ReactionID, map.reactionName(i)))*3600);
        %retrieve max presto kcat 
        enz_idx=find(raw_batchmod.S(enzMetIdx,ismember(raw_batchmod.rxns, map.final_match(i))));
        if length(enz_idx)~=1
            if length(unique(raw_batchmod.S(enzMetIdx(enz_idx),ismember(raw_batchmod.rxns, map.final_match(i)))))==1
            %complexes should not be found but accept them if kcat is
            %identical
    
            davidi_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
            
            else
                warning(['reaction with missing enzyme or enzyme complex detected comparison impossible. skipping reaction ', ...
                raw_batchmod.rxns{ismember(raw_batchmod.rxns, map.final_match(i))}])
            end
    
        else
        davidi_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
        
        end
    end
    
    ecModel_batch_keff = davidi_mod;
    
    % Run simulations for adjusted model

    fprintf('\n' + "Running simulations for  " + currentSample + " keff-adjusted model"  + '\n');
    
    tempModel = ecModel_batch_keff;
    tempSol = optimizeCbModel(tempModel, 'min', 'one');

    glucoseID = find(~cellfun('isempty',strfind(tempModel.rxns,'EX_glc__D_e_REV'))); 

    tempModel = changeRxnBounds(tempModel, 'EX_glc__D_e_REV', tempSol.x(glucoseID)*(1+0.05), 'u');
    tempModel = changeRxnBounds(tempModel, 'EX_glc__D_e_REV', tempSol.x(glucoseID)*(1-0.05), 'l');

    tempModel = changeObjective(tempModel, 'BIOMASS_Ec_iML1515_core_75p37M', 1);

    tempSol = optimizeCbModel(tempModel, 'max');

    ptotSTR = 0.61;
    f = 0.5;
    sigma = 0.5;

    etotSTR = ptotSTR*f*sigma;

    [~,nRxns] = size(ecModel_batch_keff.S);
    enzymeIds = find(~cellfun('isempty',strfind(ecModel_batch_keff.rxnNames,'prot_'))); 
    enzymeIds(end,:) = [];
    [enzymeNum,~] = size(enzymeIds);

    ecModel_batch_keff = changeRxnBounds(ecModel_batch_keff, 'BIOMASS_Ec_iML1515_core_75p37M', tempSol.f*(1+0.05), 'u');
    ecModel_batch_keff = changeRxnBounds(ecModel_batch_keff, 'BIOMASS_Ec_iML1515_core_75p37M', tempSol.f*(1-0.05), 'l');
    
    LPproblem = buildLPproblemFromModel(ecModel_batch_keff);
    
    enzymes = ecModel_batch_keff.enzymes(:);
    
    pIdx = [];
    rxnIdx = {};
    
    for i=1:numel(enzymes)
        pIdx(i,1) = find(strcmpi(ecModel_batch_keff.metNames,join(['prot_' char(enzymes(i))],"")));
        rxnIdx{i,1} = find(ecModel_batch_keff.S(pIdx(i),:) < 0);
    end
    
    % LPproblem.A(pIdx,enzymeIds) = LPproblem.A(pIdx,enzymeIds).*(etotSTR);
    % LPproblem.ub(enzymeIds) = ecModel_batch_keff.MWs;
    
    LPproblem.c = zeros(size(LPproblem.A,2),1);
    LPproblem.c(6085) = etotSTR;
    LPproblem.c(enzymeIds) = -1;
    LPproblem.c(enzymeIds) = ecModel_batch_keff.MWs.*LPproblem.c(enzymeIds);
    
    LPproblem.osense = 1;

    MinimizedFlux = solveCobraLP(LPproblem);

    MinimizedFlux.x = MinimizedFlux.full;

    % fprintf('\n');
    % fprintf('LP SOLUTION\n');
    % printFluxes(ecModel_batch_keff, MinimizedFlux.x, true);

    % Compare predicted enzyme usage to experimental proteomics data

    protModel = {};
    protModel(:,1) = ecModel_batch_keff.enzGenes;
    protModel(:,2) = ecModel_batch_keff.enzymes;
    protModel = cell2table(protModel);
    protModel.Properties.VariableNames = {'Gene' 'Protein'};
    protModel.Gene = char(protModel.Gene);
    protModel.Protein = char(protModel.Protein);

    protExp = kapp;
    toRemove = {'ReactionID', 'geneName', 'kappmax', 'subsystem_accordingToModel_'};
    protExp = removevars(protExp, toRemove);
    protExp.Properties.VariableNames = {'Gene' 'Abundance'};
    protExp.Gene = char(protExp.Gene);

    protExp = innerjoin(protExp, protModel);

    protPred = {};
    protPred(:,1) = ecModel_batch_keff.enzymes;
    protPred(:,1) = replace(protPred(:,1), 'draw_prot_', '');
    protPred(:,2) = num2cell(MinimizedFlux.x(enzymeIds));
    protPred = cell2table(protPred);
    protPred.Properties.VariableNames = {'Protein' 'Predicted'};
    protPred.Protein = char(protPred.Protein);

    protMerged = innerjoin(protExp, protPred);
    protMerged(ismember(protMerged.Predicted, 0),:)=[];
    protMerged.Abundance = protMerged.Abundance./1e12;

    protMergedLog = protMerged;
    protMergedLog.Predicted = log10(abs(protMergedLog.Predicted));
    protMergedLog.Abundance = log10(abs(protMergedLog.Abundance));

    cor = {};
    cor(1,:) = {'Pearson', 'Spearman'};
    
    [cor{2,1}, cor{3,1}] = corr(protMergedLog.Abundance, protMergedLog.Predicted, 'Type', char(cor(1,1)));
    [cor{2,2}, cor{3,2}] = corr(protMergedLog.Abundance, protMergedLog.Predicted, 'Type', char(cor(1,2)));
    
    rmse = protMergedLog;
    
    rmse.RMSE = rmse.Abundance - rmse.Predicted;
    rmse.RMSE = (rmse.RMSE).^2;
    rmseResult = median(rmse.RMSE);
    rmseResult = sqrt(rmseResult);
    
    [numEs,~] = size(protMerged);
    
    for i=1:size(protPred.Protein(:,1))
        enzMin(i,1) = convertCharsToStrings(protPred.Protein(i,:));
    end
    
    enzMinIdx = find(contains(ecModel_batch_keff.enzymes, enzMin));
    enzGenesMin = ecModel_batch_keff.enzGenes(enzMinIdx);
    % enzGenesMin = convertCharsToStrings(enzGenesMin);
    
    ecModel_batch_keffb = ravenCobraWrapper(ecModel_batch_keff);
    ecModel_batch_keffb.grRules = ecModel_batch_keffb.rules;
    ecModel_batch_keffb.geneNames = ecModel_batch_keffb.genes;
    
    rxnsEnz = findRxnsFromGenes(ecModel_batch_keffb, enzGenesMin);
    
    fieldNamesRxns = fieldnames(rxnsEnz);
    
    for i=1:numel(fieldNamesRxns)
        for j=1:size(rxnsEnz.(fieldNamesRxns{i}))
            RxnNames(i,j) = rxnsEnz.(fieldNamesRxns{i})(j);
        end
    end
    
    for m=1:size(RxnNames(:,1))
        RxnsArray.(fieldNamesRxns{m}) = RxnNames(m,:);
        keep = any(~cellfun('isempty',RxnsArray.(fieldNamesRxns{m})), 1);
        RxnsArray.(fieldNamesRxns{m}) = RxnsArray.(fieldNamesRxns{m})(:,keep);
        RxnsArray.(fieldNamesRxns{m})(:,end) = [];
    end
    
    for y=1:size(RxnNames(:,1))
        RxnsIDArray.(fieldNamesRxns{y}) = find(contains(ecModel_batch_keff.rxns, RxnsArray.(fieldNamesRxns{y})(:,:)));
    end
    
    for z=1:size(RxnNames(:,1))
        RxnsFluxArray.(fieldNamesRxns{z}) = MinimizedFlux.x(RxnsIDArray.(fieldNamesRxns{z})(:,:));
        RxnsFluxArray.(fieldNamesRxns{z}) = nonzeros(RxnsFluxArray.(fieldNamesRxns{z}));
        RxnsFluxArray.(fieldNamesRxns{z}) = max(RxnsFluxArray.(fieldNamesRxns{z}));
    end
    
    RxnsFluxArray = struct2cell(RxnsFluxArray);
    ProtMin = ecModel_batch_keff.enzymes(enzMinIdx);
    
    resultsGEM.Samples(k) = sampleIDs(k);
    resultsGEM.Pearson(k) = cor{2,1};
    resultsGEM.pvalue_P(k) = cor{3,1};
    resultsGEM.Spearman(k) = cor{2,2};
    resultsGEM.pvalue_S(k) = cor{3,2};
    resultsGEM.RMdSE(k) = rmseResult;
    resultsGEM.NumEs(k) = numEs;
    
    RxnsFlux = cell2table(cell(0,0));
    RxnsFlux.Genes = char(ecModel_batch_keffb.geneNames(enzMinIdx));
    RxnsFlux.Protein = ProtMin;
    RxnsFlux.Protein = char(RxnsFlux.Protein);
    RxnsFlux.Fluxes = RxnsFluxArray;
    RxnsFluxMerged = innerjoin(RxnsFlux, protPred);
    RxnsFluxMerged = outerjoin(RxnsFluxMerged, protExp);
    toRemove2 = {'Gene', 'Protein_protExp'};
    RxnsFluxMerged = removevars(RxnsFluxMerged, toRemove2);

    RxnsFluxFilename = "Fluxes_" + currentSample + ".xlsx";
    writetable(RxnsFluxMerged, RxnsFluxFilename)

    fprintf('\n' + "Finished for sample " + currentSample + '\n');

end

%%
resultsGEMFilename = "resultsGEM.xlsx";
writetable(resultsGEM, resultsGEMFilename)
toc
