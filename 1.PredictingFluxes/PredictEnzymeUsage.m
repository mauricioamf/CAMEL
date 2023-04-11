function MinimizedFlux = PredictEnzymeUsage(modelSTR, ptotSTR, f, sigma)
% Inputs:
%   modelSTR - a COBRA model with protein pools in the metabolites
%              (required)
%   ptotSTR - the total protein concentration (g protein / g cell)
%             (default: 0.5)
%   f - mass fraction of all measured proteins included in the model
%       (default: 0.5)
%   sigma - average in vivo enzyme saturation (default: 0.5)
%
% Output:
%   MinimizedFlux - a structure containing the solution of the LP problem

    if nargin < 4
        sigma = 0.5;
        if nargin <3
            f = 0.5;
            if nargin <2
                ptotSTR = 0.5;
            end
        end
    end

    %% Construct the LP problem to minimize ES1
    enzymeIds = find(~cellfun('isempty',strfind(modelSTR.rxnNames,'prot_')));
    enzymeIds(end,:) = [];

    LPproblem = buildLPproblemFromModel(modelSTR);

    enzymes = modelSTR.enzymes(:);

    pIdx = [];
    rxnIdx = {};

    for i=1:numel(enzymes)
        pIdx(i,1) = find(strcmpi(modelSTR.metNames,join(['prot_' char(enzymes(i))],"")));
        rxnIdx{i,1} = find(modelSTR.S(pIdx(i),:) < 0);
    end

    LPproblem.A(pIdx,enzymeIds) = LPproblem.A(pIdx,enzymeIds).*(ptotSTR*f*sigma);
    LPproblem.ub(enzymeIds) = modelSTR.MWs;

    LPproblem.c = zeros(size(LPproblem.A,2),1);
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
    MinimizedFlux.x = MinimizedFlux.full;
    MinimizedFlux.x(enzymeIds) = MinimizedFlux.x(enzymeIds);

    fprintf('\n');
    fprintf('LP SOLUTION\n');
    printFluxes(modelSTR, MinimizedFlux.x, true);

    fprintf('\n');
    fprintf('LP SOLUTION FOR ES \n');
    protNames.ids = find(~cellfun('isempty',strfind(modelSTR.metNames,'prot_pool')));
    protNames.names = modelSTR.metNames(protNames.ids);
    printFluxes(modelSTR, MinimizedFlux.x, false, '', '', '', protNames.names);
end
