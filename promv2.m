function [f, f_KO, v, v_KO, status1, lostXns, probTFgene] =  promv2...
          (model, expression, expressionID, regulator, targets, ...
          litEvidence, probPrior, subsets, minFlux, maxFlux, KAPPA, ...
          dataThreshVal, probTFgene, sizeFlag)

% The PROM algorithm predicts the growth phenotype and the flux response
% after transcriptional perturbation, given a metabolic and regulatory
% network.
% 
% USAGE:
%   
%    [f, f_KO, v, v_KO, status1, lostXns, probTFgene] =  promv2...
%    (model, expression, expressionID, regulator, targets, ...
%    litEvidence, probPrior, subsets, minFlux, maxFlux, KAPPA, ...
%    dataThreshVal, probTFgene, sizeFlag)
%
% INPUTS:
%    model:           The metabolic model obtained from COBRA toolbox 
%                     through the readcbmodel command. 
%    expression:      Gene expression data. No need to normalize or impute.
%                     rows = genes, columns = conditions 
%    expressionID:    An array of identifiers for each row/gene.
%    regulator:       Regulatory network. Format includes a cell array of 
%                     regulators and matching target genes.
%                     Example:
%                         regulator = {'RegA'; 'RegB' ; 'RegC'};  
%                         targets ={'GeneA'; 'GeneB'; 'GeneC'}
%                     Note that the names or identifiers used in the 
%                     regulatory data should match the names/IDs given for 
%                     gene expression data.
%    minFlux:         Obtained either from fastFVA or fluxvariability 
%                     command in COBRA.
%    maxFlux:         Obtained either from fastFVA or fluxvariability 
%                     command in COBRA.
%                     Example: [minFlux, maxFLux] = fastFVA(model);
%    sizeFlag:        Tells PROM if the regulatory network is large. It is  
%                     0 for large networks and 1 for small networks (less 
%                     than 1000 interactions).
%
% OPTIONAL INPUTS:
%    litEvidence:     High confidence interactions (not necessarily based  
%                     on literature) should be flagged as 1 in litEvidence
%                     array and the rest should be set as 0. Should have 
%                     same length as regulator/target arrays.  
%    probPrior:       Should be set to values between 0 and 1 for those 
%                     interactions with literature evidence (other values 
%                     in the array would be ignored). Should have same 
%                     length as regulator/target arrays.    
%    KAPPA:           Determines strength of regulatory constraints.  
%                     Default value of 1 works for most systems.
%    subsets:         Subsets of transcription factors (TFs) for which prom
%                     should be run. Default includes all TFs.
%
% OUTPUTS:
%    f:               Growth rate. 
%    f_KO:            Growth rate after regulator knockout.
%    v:               Flux response.
%    v_KO:            Flux response after regulator knockout.
%    status1:         The glpk solver status. Should be 5 for glpk (if it's
%                     not, check solver error log.
%    lostXns:         Gives the interactions that could not be quantified 
%                     based on the threshold set for binarization. The 
%                     program would shoot a warning if the threshold is 
%                     poorly chosen. The default value (0.2 - 0.4) should 
%                     work for most cases.
%    probTFgene:      Gives the probabilities estimated for each 
%                     interaction.
%
% EXAMPLE:
%
%    load mtbpromdata    % tuberculosis model used in PNAS paper
%    f_KO = promv2(model, expression, expressionid, regulator, targets);
%    [minFlux, maxFlux] = fastFVA(model);
%    [f, f_KO, v, v_KO, status1, lostXns, probTFgene] =  promv2(model, ...
%     expression, expressionid, regulator, targets, [], [], [], ...
%    minFlux, maxFlux, [], [], []);
%
% NOTE:
%    ***This is a very important information to be highlighted***
%
% .. Author: - Sriram Chandrasekaran, 11/1/20, University of Michigan
%%
if nargin == 5     % check for correct number of inputs
     litEvidence = [];
     probPrior = [];
     subsets = [];
elseif (nargin < 5) || (nargin == 6)
    error('Incorrect number of input arguments to %s',mfilename)
end

% Variable initialization 
disp('initializing data')
regulated = targets;
TFnames = unique(regulator);    % transcription factor names    
weights = model.c;   % linear objective
stoich  = model.S;   % stoichiometric matrix
cType = repmat('S',size(model.b));
lbFF = model.lb;    % reaction lower bounds
ubFF = model.ub;    % reaction upper bounds
dxdt = model.b;     % accumulation of metabolites
param.msglev = 1;
lbFF(lbFF == ubFF) = ubFF(lbFF == ubFF) - 1E-6;

if (~exist('KAPPA','var')) || (isempty(KAPPA))         
    KAPPA = 1;    % assign 1 to KAPPA if no input provided
end

% Flux Variability Analysis

% for small networks its faster to compute min/max for the subset 
% that's regulated rather than the entire metabolic model
if (~exist('sizeflag', 'var')) || (isempty(sizeFlag))
    if (length(regulator) < 1000) 
        sizeFlag = 1;
    else
        sizeFlag = 0;
    end
end

if ~sizeFlag
    if (~exist('minFlux', 'var')) || (isempty(minFlux))  
        [minFlux, maxFlux] = fluxvariability(model); 
        %[minFlux, maxFlux] = fastFVA(model); 
    end
end

if (~exist('dataThreshVal','var')) || (isempty(dataThreshVal))
    dataThreshVal = 0.33; 
end


[rxnPos, geneList] = find(model.rxnGeneMat); % find rxn positions 
BnumsInExpsn = expressionID;
litEvidence = logical(litEvidence);

% BnumsToBeKO gives the geneIDs of the genes to be knocked out one by one
% if subset of TFs not provided, default to all TFs
if ~isempty(subsets)
    BnumsToBeKO = subsets;
    regulated = targets(ismember(regulator, subsets));
    regulator = regulator(ismember(regulator, subsets));
else
    BnumsToBeKO= TFnames;   
end

scou = 1;
%% Additional variable initializations
lbG = model.lb; 
ubG = model.ub;
lbG(lbG == ubG) = ubG(lbG == ubG) - 1E-6;
a1 = [stoich, zeros(size(stoich, 1), length(ubG)), zeros(size(stoich, 1),...
     length(ubG))];
a2 = sparse([eye(length(ubG)), eye(length(ubG)), zeros(length(ubG))]);
a3 = sparse([eye(length(ubG)), zeros(length(ubG)), -eye(length(ubG))]);
A = [a1; a2; a3];
weights11 = [weights; zeros(2 * length(lbG), 1)];
weights00 = [weights; zeros(2 * length(lbG), 1)];

lb11 = [-1000 * ones(length(lbG), 1);zeros(length(lbG), 1);...
    zeros(length(lbG), 1)];
ub11 = [1000 * ones(length(lbG),1);zeros(length(lbG),1);...
    zeros(length(lbG), 1)];
lb11(lb11 == ub11) = ub11(lb11 == ub11) - 1E-6;

dxdt0 = [zeros(size(stoich, 1),1); lbG; ubG];
ctype1 = [repmat('S', size(stoich, 1),1); repmat('L', size(lbG, 1), 1);...
    repmat('U', size(lbG, 1), 1)];
[v0, f0] = glpk(-weights11, A, dxdt0, lb11, ub11, ctype1, [], [], param);

%% Find Probabilities using a Global Threshold
tic
lost_xn = false(size(regulated));
remove_interactions1 = false(size(regulated));
disp('finding probabilities')
counter = 1; 
counter1 = 1; 
expressionProcessed = expression;
expressionProcessed = knnimpute(expressionProcessed);
expressionProcessed = quantilenorm(expressionProcessed);  %its already normalized..
expressionNoThresh = expressionProcessed;

dataThresh = quantile(expressionProcessed(:), dataThreshVal);
if dataThresh < 0
    expressionProcessed(expressionProcessed >= dataThresh) = 1;
    expressionProcessed(expressionProcessed < dataThresh) = 0;
else
    expressionProcessed(expressionProcessed < dataThresh) = 0;
    expressionProcessed(expressionProcessed >= dataThresh) = 1;
end

if (~exist('probTFgene','var')) || (isempty(probTFgene))
    try probTFgene = ones(size(regulated));
    catch M1
        probTFgene = 1; 
    end
    for i = 1:length(regulated)
        k = find(ismember(BnumsInExpsn, regulated(i)));
        l = find(ismember(BnumsInExpsn, regulator(i)));
        if ~isempty(k) && ~isempty(l)
            te = expressionNoThresh(k,:);
            te1 = expressionNoThresh(l,:);
            tec = expressionProcessed(k,:);
            tec1 = expressionProcessed(l,:);
            counter1 = counter1 + 1;
            try kstest2(te(tec1 == 1), te(tec1 == 0));
                if  (kstest2(te(tec1 == 1), te(tec1 == 0)) == 1)
                    prob1 = sum(tec(tec1 == 0)) / length(tec(tec1 == 0));
                    probTFgene(i) = prob1;
                    counter = counter + 1;
                else
                    probTFgene(i) = 1;    % no effect
                end
            catch ERRLG
                % can't be estimated from microarray (if it has strong 
                % evidence, consider setting this to zero)
                probTFgene(i) = 1;
                lost_xn(i) = 1;
            end
        else
            probTFgene(i) = 1;
            lost_xn(i) = 1;
        end
    end
    probTFgene = probTFgene(:);    % flatten array
    toc
    % Uncomment below condition to check for issue with binarization (if
    % over 75% interactions missed, should change threshold)
%    if (sum(lost_xn) > 0.75*length(probtfgene))  
%        % missed then its time to change the threshold
%        datathreshflag = 1;
%        disp('change binarization threshold')
%    else
    dataThreshFlag = 0;
%    end
    
%   you could set those interactions that you think have strong literature 
%   evidence to have predefined probabilities
    if ~isempty(litEvidence)
        probTFgene(litEvidence) = probPrior(litEvidence);  
    end                                                            
    toc
else
    disp('Using pre-supplied probabilities')
    dataThreshFlag = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run PROM for each knockout
disp('Running PROM')
thresh = 10^(-6); 
mThresh = 10^(-3);
allGenes = [model.genes; TFnames];
[ir, posGeneList] = ismember(regulated, model.genes);
lbf = lbFF;
ubf = ubFF;
[v, f_wt] = glpk(-weights, stoich, dxdt, lbf, ubf, cType);
f_wt = -f_wt;

v(abs(v) < thresh) = 0; % floor small values to zero
kappa = KAPPA;
hw = waitbar(0);

% knock out the genes
vm = zeros(size(v));
if ~sizeFlag
    % disp('running fast fva for the entire network'); tic;
    %         [v11,v12] = fastFVA(model); toc;
    minFlux(abs(minFlux) < thresh) = 0;
    maxFlux(abs(maxFlux) < thresh) = 0;
end

for geneCount = 1:length(BnumsToBeKO) 
    disp(geneCount)
    lbG = lbf; 
    ubG = ubf;
    lb11 = [-1000 * ones(length(lbG), 1); zeros(length(lbG), 1);...
        zeros(length(lbG), 1)];
    ub11 = [1000 * ones(length(lbG), 1); zeros(length(lbG), 1);...
        zeros(length(lbG), 1)];
    
    % check if gene is metabolic, regulatory or both
    if any(strcmpi(model.genes, BnumsToBeKO(geneCount)))
        tempPos = rxnPos(geneList == find(strcmp(model.genes,...
            BnumsToBeKO(geneCount))));
        for jj = 1:length(tempPos)
            if model.rev(tempPos(jj))
                lbG(tempPos) = -thresh;
                ubG(tempPos) = thresh;
            else
                lbG(tempPos) = -thresh;
            end
        end
    end
    
    [v1, fk(geneCount)]  = glpk(-weights, stoich, dxdt, lbG, ubG, cType);
    v1(abs(v1) < thresh) = 0;   % floor small values to zero

    if any(ismember(TFnames, BnumsToBeKO(geneCount)))
        tfState = logical(zeros(size(TFnames)));
        tfState(find(ismember(TFnames, BnumsToBeKO(geneCount)))) = 1;
        k = find(ismember(regulator, TFnames(tfState)));
        tempGene = regulated(k);
        tempGeneProbs = probTFgene(k);
        tempGenePos = posGeneList(k);
        tempRxnPos = rxnPos(ismember(geneList, tempGenePos));

        % gene-protein-reaction relationship
        x = true(size(model.genes));
        [isInModel, geneInd] = ismember(tempGene, model.genes);
        x(geneInd) = false;
        constrainRxn = false(length(tempRxnPos),1);
        % figure out if any of the reaction states have changed
        for j = 1:length(tempRxnPos)
            if (~eval(model.rules{tempRxnPos(j)}))
                constrainRxn(j) = true;
            end
        end
        
        % Constrain flux through the reactions associated with these genes
        tempGeneProbs(tempGenePos == 0)  = '';
        tempGene(tempGenePos == 0) = '';
        tempGenePos(tempGenePos == 0)  = '';
        
        % tempRxnPos has the rxns that are going to be affected by this TF.
        % kRxnPos are the rxns that will be affected by this target gene alone.
        % Loop around all the genes.
        for l = 1:length(tempGenePos)
            if ~isnan(tempGeneProbs(l))
                kRxnPos = ismember(tempRxnPos, rxnPos(ismember(geneList,...
                    tempGenePos(l))));
                for m = 1:length(tempRxnPos)
                    if kRxnPos(m)
                        if constrainRxn(m)
                            % if tempGeneProbs is 1, no use in changing the
                            % bounds 
                            if (tempGeneProbs(l) < 1)
                                % if it's zero, no point in estimating vm 
                                % again, but cannot include in the above 
                                %statement because must change the bounds
                                if (tempGeneProbs(l) ~= 0)   
                                    % done to save time -> if estimated \
                                    % already, use it
                                    if ~vm(tempRxnPos(m))    
                                        if sizeFlag
                                            weights1 = weights; 
                                            lbv = lbf; 
                                            ubv = ubf;
                                            grwthpos = find(weights == 1);
                                            lbv(grwthpos) = v(grwthpos);
                                            weights1(tempRxnPos(m)) = -1;
                                            [minFlux, fva1]  = glpk(-weights1, stoich, dxdt, lbv, ubv, cType);
                                            weights1(tempRxnPos(m)) = 1;
                                            [maxFlux, fva2]  = glpk(-weights1, stoich, dxdt, lbv, ubv, cType);
                                        end
                                        if v(tempRxnPos(m)) < 0
                                            vm(tempRxnPos(m)) = min(...
                                                [minFlux(tempRxnPos(m)),...
                                                maxFlux(tempRxnPos(m)),...
                                                v(tempRxnPos(m))]);
                                        elseif v(tempRxnPos(m)) > 0
                                            vm(tempRxnPos(m)) = max(...
                                                [minFlux(tempRxnPos(m)),...
                                                maxFlux(tempRxnPos(m)),...
                                                v(tempRxnPos(m))]);
                                        else
                                            vm(tempRxnPos(m)) = max(...
                                                [abs(minFlux(tempRxnPos(m))),...
                                                abs(maxFlux(tempRxnPos(m))),...
                                                abs(v(tempRxnPos(m)))]);
                                        end
                                    end
                                    %end
                                end
                                % flux times probability
                                fluxProb = vm(tempRxnPos(m)) * tempGeneProbs(l); 
                                % Make sure we aren't violating the 
                                % original bounds.
                                % Also get the lowest value if there were 
                                % multiple modifications for the rxn.
                                if  v(tempRxnPos(m)) < 0 
                                    tem = max([lbf(tempRxnPos(m)), fluxProb,...
                                        lbG(tempRxnPos(m))]); 
                                    % prevents the solver from crashing
                                    lbG(tempRxnPos(m)) = min([tem, -thresh]);    
                                    ub11(1 * length(ubG) + tempRxnPos(m))...
                                        = 1000;    % unconstrain model
                                    
                                    weights11(1 * length(ubG) + tempRxnPos(m))...
                                        = (-1 * kappa / abs(vm(tempRxnPos(m))))...
                                        * abs(f0);   % v0, f0 are the wild type values..
                                    vv = max([abs(vm(tempRxnPos(m))), mThresh]);

                                    weights11(1 * length(ubG) + ...
                                        tempRxnPos(m)) = min([(kappa * (-1)...
                                        * abs(f0)) / (abs(vv)),...
                                        weights11(1*length(ubG) + tempRxnPos(m))]);
                                elseif v(tempRxnPos(m)) > 0
                                    tem = min([fluxProb, ubf(tempRxnPos(m)), ubG(tempRxnPos(m))]);
                                    ubG(tempRxnPos(m)) = max(tem, thresh);
                                    ub11(2 * length(ubG) + tempRxnPos(m)) = 1000;
                                    vv = max([abs(vm(tempRxnPos(m))), mThresh]);
                                    % new weights based on kappa, 
                                    % normalized with growth rate
                                    weights11(2 * length(ubG) + tempRxnPos(m))...
                                        = min([(kappa * (-1) * abs(f0))...
                                        / abs(vv), weights11(2 * length(ubG)...
                                        + tempRxnPos(m))]);  
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    dxdt0 = [zeros(size(stoich, 1), 1); lbG; ubG];
    param.itlim = 1000000;
    
    % optimizeCbModel
    [v00(geneCount, :), f00(geneCount), status1(geneCount)] = glpk(...
        -weights11, A, dxdt0, lb11, ub11, ctype1, [], [], param);
    f(geneCount) = v00(geneCount, find(weights));
    counter  = 1;
    lbh = lbG;
    ubh  = ubG;
    
    while ((status1(geneCount) ~= 5) || (v00(geneCount, find(weights)) < 0))
        lbh(lbh ~= lbf) = lbh(lbh ~= lbf) - 1E-3;
        ubh(ubh ~= ubf) = ubh(ubh ~= ubf) + 1E-3;
        dxdt0 = [zeros(size(stoich, 1), 1); lbh; ubh];
        [v00(geneCount, :), f00(geneCount), status1(geneCount)] = glpk(...
            -weights11, A, dxdt0, lb11, ub11, ctype1, [], [], param);
        counter = counter + 1
        if (counter > 50)
            % if it's impossible to estimate, then
            % check the unweighted g.r and the one with max weight prom -
            % if much less difference use it - any way warn the user about
            % the problem at the iteration number;
            [v3 ,f3, status3] = glpk(-weights, stoich, dxdt, lbG, ubG, cType);
            [v30, f30, status30] = glpk(-weights00, A, dxdt0, lb11, ub11, ctype1);
            %if abs((f3- v30(find(weights)))/abs(f3)) < 0.1
            f00(geneCount) = -f3;
            v00(geneCount, find(weights)) = -f3;
            %else
            disp(' problem in'); 
            disp(geneCount);
            break;  
            %end
            disp('check'); 
            disp(geneCount); 
            break;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f(geneCount) = v00(geneCount, find(weights));

    lbg_st(geneCount, :) = lbG;
    ubg_st(geneCount, :) = ubG;
    
    lb_st(geneCount, :) = lbFF;
    ub_st(geneCount, :) = ubFF;
    
    [v2(geneCount, :), f1(geneCount), status] = glpk(-weights, stoich,...
        dxdt, lbG, ubG, cType);
    f_KO(geneCount) = -f1(geneCount);
    
    kTime = toc;
    waitPar = [num2str(ceil(geneCount / length(TFnames) * 100)),...
        '% complete. Time taken:', num2str(ceil(kTime)),' secs'];
    %waitbar(ci/length(tfnames),hw,waitpar);
    
    ff00(scou, geneCount) = v00(geneCount,find(weights));
    %disp(ff00(scou,ci))
    
    if dataThreshFlag
        % if none of the probabilities of a gene can be estimated:
        if all(lost_xn(k))    
            disp('Error: xns cant be estimated')
            v00(geneCount,:) = NaN;
            f1(geneCount) = NaN;
            v2(geneCount,:) = NaN;
            f00(geneCount) = NaN;
            %break;
        end
    end
    
    
    clear tempGenePos tempGeneProbs tempRxnPos k
    %disp(bnumstobekoed(ci))
end

% f_ko = -f1';
v_KO = v2;

% f00_ko(:,scou) = v00(:,find(weights));
v00_ko = v00;

% f = f00_ko;
v = v00_ko;

lostXns(:,scou) = lost_xn;

%disp([f,f_ko])

end

