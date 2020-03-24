function [f,initial_network,final_network] =  GEMINI(model,expression,expressionid,regulator,targets,phenotype,type,subsets,v11,v12,sizeflag,OPTIMAL_THRESH,metric_type,DATATHRESHVAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [f,initial_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype,type,subsets,v11,v12);
% GEMINI performs iterative refinement of regulatory networks using
% metabolic networks
% INPUTS
% Model is the metabolic model for the organism (obtained from COBRA
% toolbox through readcbmodel command) . The model should be set to a
% specific growth condition under study ( like glucose minimal media)
%
% Gene expression data - rows - genes,columns - conditions; (preferably
% normalized and imputed)
%
% Expressionid - an array of identifiers for each row/gene should be included
%
% draft regulatory network - format - cell array of regulators and matching target genes
% example   Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets =
% {'GeneA';'GeneB';'GeneC'}
% note that the names or identifiers used in the regulatory data should
% match the names/ids given for gene expression data
%
% phenotype - logical vector (true/false) - the growth phenotype of each transcription factor knockout under a specific condition specified
% by the metabolic model. it should be the same length as the number of
% TFs in the model.  true (1) - lethal; false (0) - viable
%
% type - a string describing the knockout phenotype data - should be either 'lethal' (default) or 'suboptimal' ;
% this is required to set the OPTIMAL_THRESH value; OPTIMAL_THRESH is the threshold for determining lethal/non-lethal phenotypes; by default, a knockout
% that grows less than 5% of the wild type growth rate is considered lethal
% and less than 95% of the wildtype growth rate is considered suboptimal
% ; ( OPTIMAL_THRESH default value = 0.05 for 'lethal' type, corresponds to 5% and 0.95 for 'suboptimal')
%
%OPTIONAL
% % subsets : subsets of tfs for which GEMINI should be run ; default - for all tfs
%v11,v12 are the minimum and maximum possible flux through each reaction in the model for the given
% condition. this obtained through flux variability analysis (either from fastfva or fluxvariability command in
%cobra)
% >> [v11,v12] = fluxvariability(model);
%  sizeflag tells GEMINI if the regulatory network is large. it is 0 for large networks
%  and 1 for small networks ( less than 1000 interactions)
%
% OUTPUT - the algorithm gives the refined network and the growth rate (f) after knock out of all
% tfs in the refined network.
% note that f is only semi-quantitative, unless trained on suboptimal growth data
% EXAMPLE
% load yeast_gemini_data regulator targets expression expressionid model v11 v12 phenotype
% [f,initial_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype,'lethal',{'YAL051W'},v11,v12);

%===========================================================
%% INPUT HANDLING
%===========================================================
if (~exist('sizeflag','var')) || (isempty(sizeflag))
    if (length(regulator) < 1000) % for small networks its faster to compute min/max for the subset thats regulated than the entire metab model
        sizeflag = 1;
    else
        sizeflag = 0;
    end
end

if ~sizeflag
    if (~exist('v11','var')) || (isempty(v11))
        [v11,v12] = fluxvariability(model); % or use fastfva
        %[v11,v12] = fastFVA(model);
    end
end

if (~exist('OPTIMAL_THRESH','var')) || (isempty(OPTIMAL_THRESH))
if (~exist('type','var')) || (isempty(type))
    OPTIMAL_THRESH = 0.05;
    disp('phenotype data assumed to be lethal/non lethal:optimal_thresh set to 0.05')
else
    type = upper(type); type = strtrim(type);
    if ismember(type,{'LETHAL'})
        OPTIMAL_THRESH = 0.05;
    elseif ismember(type,{'SUBOPTIMAL'})
        OPTIMAL_THRESH = 0.95;
    else
        error('incorrect input for phenotype data type. should be either "lethal" or "suboptimal"')
    end
end
end


if (~exist('DATATHRESHVAL','var') || isempty(DATATHRESHVAL))
     DATATHRESHVAL = 0.33;
end

if (~exist('metric_type','var')|| isempty(metric_type))

    metric_type = 'default';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%kappavec = [0.01,0.1,1,5,10,100,1000,10000];
kappavec  = 10;
for iter = 1:length(kappavec),
    % for each iteration set kappa for tuning
    
    KAPPA = kappavec(iter);
    
    fprintf('params used - KAPPA: %d and OPTIMAL_THRESH: %1.3f \n', KAPPA,OPTIMAL_THRESH)
    
    [f(iter,:),remov_interactions] =  GEMINI_phenotype(model,expression,expressionid,regulator,targets,phenotype,subsets,v11,v12,OPTIMAL_THRESH,KAPPA,sizeflag,metric_type,DATATHRESHVAL);
    initial_network = strcat(regulator,'--',targets);
    final_network(:,iter)  = initial_network(~remov_interactions); %#ok<AGROW>
    
    disp('number of interactions removed for this TF:')
    disp(sum(remov_interactions))
end

end
