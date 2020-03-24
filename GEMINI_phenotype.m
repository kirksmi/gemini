function [f,remove_interactions1] =  GEMINI_phenotype(model,expression,expressionid,regulator,targets,phenotype,subsets,v11,v12,OPTIMAL_THRESH,KAPPA,sizeflag,metric_type,DATATHRESHVAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [f,f_ko,v,v_ko,status1,lostxns,remove_interactions1] =
%% GEMINI_phenotype(model,expression,expressionid,regulator,targets,phenotype,litevidence,prob_prior,subsets,v11,v12,OPTIMAL_THRESH,KAPPA,datathresh)
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
% by the metabolic model
%
% KAPPA - determines strength of regulatory constraints ; this value is
% data-specific and obtained by tuning over a range of values from 0.1-1000


%OPTIONAL
% subsets - subsets of tfs for which prom should be run ; default - run for all tfs

% OUTPUT - the algorithm gives the growth rate (f) and flux response (v) after knock out of all
% the regulators in the regulatory model; status is the glpk solver status
% the status should be 5 for glpk; if its not then check solver error log
% lostxns gives the interactions that could not be quantified based on the
% threshold set for binarization
% the program would shoot a warning if the threshold chosen is bad. The
% default value (0.2 - 0.4) should work for most cases
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INPUT HANDLING
%===========================================================
%SOME BASIC INITIALIZATION
%===========================================================
if ~exist('metric_type','var')
    metric_type = 'default';
end


disp('initializing data')
regulated = targets;
[tfnames,b,m] = unique(regulator);
weights = model.c; stoic  = model.S; S = model.S; ctype = repmat('S',size(model.b));lbff = model.lb; ubff = model.ub;dxdt = model.b; param.msglev = 1;
lbff(lbff==ubff) = ubff(lbff == ubff) - 1E-6;


if (~exist('DATATHRESHVAL','var') || isempty(DATATHRESHVAL))

    DATATHRESHVAL = 0.33;

end

if isempty(OPTIMAL_THRESH)
    OPTIMAL_THRESH = 0.05;
end

if (length(OPTIMAL_THRESH) ~= length(tfnames)) % that means only one value was supplied
    OPTIMAL_THRESH = repmat(OPTIMAL_THRESH(1),size(tfnames));
end

if isempty(KAPPA)
    KAPPA = 1;
end


% fprintf('params used - KAPPA: %d and DATATHRESH: %1.3f \n', KAPPA,DATATHRESHVAL)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u,v] = find(model.rxnGeneMat);
% u - reaction; v - genes
% finding rxn position
rxnpos = u;
genelist = v;
clear u v

% i need to find the reactions that each gene has influence on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bnumsinexpsn = expressionid;
%litevidence = logical(litevidence);
if ~isempty(subsets)
    bnumstobekoed = subsets;
    regulated = targets(ismember(regulator,subsets));
    regulator = regulator(ismember(regulator,subsets));
    
else
    bnumstobekoed= tfnames;   % bnumstobekoed -  gives the geneids of the genes to be knocked out - by default it knocksout all the tfs in the model one by one
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scou = 1;
%% new additions

lbg = model.lb; ubg = model.ub;
lbg(lbg==ubg) = ubg(lbg == ubg) - 1E-6;
a1 = [S,zeros(size(S,1),length(ubg)),zeros(size(S,1),length(ubg))];
a2 = sparse([eye(length(ubg)), eye(length(ubg)),zeros(length(ubg))]);
a3 = sparse([eye(length(ubg)), zeros(length(ubg)),-eye(length(ubg))]);
A = [a1;a2;a3];
weights11 = [weights;zeros(2*length(lbg),1)];
weights00 = [weights;zeros(2*length(lbg),1)];
lb11 = [-1000*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];
ub11 = [1000*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];
lb11(lb11==ub11) = ub11(lb11 == ub11) - 1E-6;
dxdt0 = [zeros(size(S,1),1);lbg;ubg];
ctype1 = [repmat('S',size(S,1),1);repmat('L',size(lbg,1),1);repmat('U',size(lbg,1),1)];
[v0,f0] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
disp('optimizing model. wild type growth rate:')
disp(v0(find(weights)))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Probabilities using a Global Threshold
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
lost_xn = false(size(regulated));
remove_interactions1 = false(size(regulated));
disp('finding probabilities')
cou = 1;cou1 = 1;cou3 = 1;
%data= knnimpute(datbackup);
data = expression;
data = knnimpute(data);
data = quantilenorm(data);  %its already normalized..
data1 = data;
%if isempty(datathresh)
datathresh = quantile(data(:),DATATHRESHVAL);
%end

%datathresh = quantile(data(:),DATATHRESHVAL);
if datathresh < 0,
    data(data>=datathresh) = 1;
    data(data < datathresh) = 0;
else
    data(data < datathresh) = 0;
    data(data>=datathresh) = 1;
end

for  i = 1:length(regulated)
    k = find(ismember(bnumsinexpsn,regulated(i)));
    l = find(ismember(bnumsinexpsn,regulator(i)));
    if ~isempty(k) & ~isempty(l)
        te = data1(k,:);
        te1 = data1(l,:);
        
        tec = data(k,:);
        
        tec1 = data(l,:);
        
        cou1 = cou1 + 1;
        try kstest2(te(tec1 == 1),te(tec1== 0));
            
            if  (kstest2(te(tec1 == 1),te(tec1== 0)) == 1),
                
                prob1 = sum(tec(tec1 == 0))/length(tec(tec1==0));
                probtfgene(i) = prob1;
                cou = cou + 1;
                %   this formula also gives the same answer  - (sum(~tec1(tec == 1))/length(tec1(tec==1))) * (sum(tec)/sum(~tec1))
                
            else
                
                probtfgene(i) = 1;  % no effect
                
            end
            
        catch ERRLG    % cant be estimated from microarray ; if it has strong evidence, i might consider setting this to zero later on
            probtfgene(i) = 1;
            lost_xn(i) = 1;
            
        end
    else
        probtfgene(i) = 1;
        lost_xn(i) = 1;
    end
end

probtfgene = probtfgene(:);
toc
if (sum(lost_xn) > 0.75*length(probtfgene))   % check if there is problem with binarizartion.. usually there ll be some interactions that'd be missed, but if most of it (like >75% ) are...
    % missed then its time to change the threshold
    datathreshflag = 1;
    disp('warning: change binarization threshold')
else
    datathreshflag = 0;
end

% if ~isempty(litevidence)
%     probtfgene(litevidence) = prob_prior(litevidence);  % u could set those interactions that u think have strong literature evidence to have predefined
% end                                                            % probabilities

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  run PROM for each knockout
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('running PROM')
count = 1; thresh = 10^(-6); mthresh = 10^(-3);
allgenes = [model.genes;tfnames];
[ir,posgenelist] = ismember(regulated,model.genes);
lbf = lbff; ubf = ubff;
[v,f_wt] = glpk(-weights,stoic,dxdt,lbf,ubf,ctype);
f_wt = -f_wt;
% flooring small values to zero
v(abs(v) < thresh) = 0;
% kapavals = [0,0.0001,0.001,0.05,0.1,0.25,0.33,0.5,1,5,10];
%kappa = kapavals(scou);
kappa = KAPPA;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hw = waitbar(0);
vm = zeros(size(v));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
weights1 = weights; lbv = lbf; ubv = ubf;
% flooring small values to zero
v11(abs(v11) < thresh) = 0;
v12(abs(v12) < thresh) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ci = 1:length(bnumstobekoed)
    
    disp(ci)
    
    lbg = lbf; ubg = ubf;
    lb11 = [-1000*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];
    ub11 = [1000*ones(length(lbg),1);zeros(length(lbg),1);zeros(length(lbg),1)];
    % check if its a metabolic or regulatory gene or both
    
    if any(strcmpi(model.genes,bnumstobekoed(ci)))
        temppos = rxnpos(genelist == find(strcmp(model.genes,bnumstobekoed(ci))));
        for jj = 1:length(temppos)
            if model.rev(temppos(jj))
                lbg(temppos) = -thresh;
                ubg(temppos) = thresh;
            else
                lbg(temppos) = -thresh;
            end
        end
        
    end
    
    
    
    [v1,fk(ci)]  = glpk(-weights,S,dxdt,lbg,ubg,ctype);
    % flooring small values to zero
    v1(abs(v1) < thresh) = 0;
    
    if any(ismember(tfnames,bnumstobekoed(ci))),
        
        
        tfstate = logical(zeros(size(tfnames)));
        tfstate(find(ismember(tfnames,bnumstobekoed(ci)))) = 1;
        k = find(ismember(regulator,tfnames(tfstate)));
        %    k(tempgeneprobs == 1) = '';
        tempgene = regulated(k);
        tempgeneprobs = probtfgene(k);
        tempgenepos = posgenelist(k);
        temprxnpos = rxnpos(ismember(genelist,tempgenepos));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % this section is for gene-protein-reaction relationship
        x = true(size(model.genes));
        [isInModel,geneInd] = ismember(tempgene,model.genes);
        x(geneInd) = false;
        constrainRxn = false(length(temprxnpos),1);
        % Figure out if any of the reaction states is changed
        for j = 1:length(temprxnpos)
            if (~eval(model.rules{temprxnpos(j)}))
                constrainRxn(j) = true;
            end
        end
        % Constrain flux through the reactions associated with these genes
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tempgeneprobs(tempgenepos == 0)  = '';
        tempgene(tempgenepos == 0) = '';
        tempgenepos(tempgenepos == 0)  = '';
        
        % temprxnpos has the rxns that are going to  be affected by this tf
        % krxnpos are the rxns that will be affected by this target gene alone..
        % we loop around all the genes..
        
        
        
        for l = 1:length(tempgenepos)
            if ~isnan(tempgeneprobs(l))
                krxnpos = ismember(temprxnpos,rxnpos(ismember(genelist,tempgenepos(l))));
                for m = 1:length(temprxnpos)
                    if krxnpos(m)
                        
                        if constrainRxn(m)
                            if (tempgeneprobs(l) < 1)    % if its 1 no use in changing the bounds - might as well save time
                                if (tempgeneprobs(l) ~= 0)   % if its zero no point in estimating vm again - saves time.. but cant include in the above statement coz u have to change the bounds
                                    %if v(temprxnpos(m))
                                    if ~vm(temprxnpos(m))    % done to save time - if estimated already use it
                                        if  v(temprxnpos(m)) < 0
                                            
                                            vm(temprxnpos(m)) = min([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        elseif v(temprxnpos(m)) > 0
                                            vm(temprxnpos(m)) = max([v11(temprxnpos(m)),v12(temprxnpos(m)),v(temprxnpos(m))]);
                                        else
                                            vm(temprxnpos(m)) = max([abs(v11(temprxnpos(m))),abs(v12(temprxnpos(m))),abs(v(temprxnpos(m)))]);
                                        end
                                        
                                        
                                    end
                                    %end
                                end
                                
                                
                                xx =    vm(temprxnpos(m))*tempgeneprobs(l); % flux times probability
                                % xx =  ( vm(temprxnpos(m))*(1 - (1/tempsubsysmid) ) + ( vm(temprxnpos(m))*tempgeneprobs(l)*(1/tempsubsysmid)));
                                if  v(temprxnpos(m)) < 0
                                    
                                    tem = max([lbf(temprxnpos(m)),xx,lbg(temprxnpos(m))]);  %make sure we arent violating the original bounds; also get the lowest value if there were multiple modifications for the rxn
                                    lbg(temprxnpos(m)) = min([tem,-thresh]);   % prevents the solver from crashing
                                    ub11(1*length(ubg) + temprxnpos(m)) = 1000;
                                    weights11(1*length(ubg) + temprxnpos(m)) = (-1*kappa/abs(vm(temprxnpos(m))))*abs(f0);   % v0 f0 are the wild type values..
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    weights11(1*length(ubg) + temprxnpos(m)) = min([(kappa*(-1)*abs(f0))/(abs(vv)), weights11(1*length(ubg) + temprxnpos(m)) ]);
                                elseif v(temprxnpos(m)) > 0
                                    tem = min([xx,ubf(temprxnpos(m)),ubg(temprxnpos(m))]);
                                    ubg(temprxnpos(m)) = max(tem,thresh);
                                    ub11(2*length(ubg) + temprxnpos(m)) = 1000;
                                    vv = max([abs(vm(temprxnpos(m))),mthresh]);
                                    weights11(2*length(ubg) + temprxnpos(m)) = min([(kappa*(-1)*abs(f0))/abs(vv), weights11(2*length(ubg) + temprxnpos(m)) ]);  % new weights based on kappa, normalized with growth rate
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    dxdt0 = [zeros(size(S,1),1);lbg;ubg];
    %lb11(lb11==ub11) = ub11(lb11 == ub11) - 1E-6; % prevents solver from crashing
    param.itlim = 1000000;
    %  optimizeCbModel
    [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
    f(ci) = v00(ci,find(weights));
    disp('predicted growth rate:')
    disp(f(ci))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% comparing phenotype to edit network
    growth(ci) = (f(ci) < OPTIMAL_THRESH(ci)*f_wt); % viable
    remove_interactions = false(size(probtfgene));
    if ~isnan(phenotype(ci))
        if (growth(ci) ~= phenotype(ci)) % if they dont match change the network
            disp('phenotype mismatch')
            if (phenotype(ci) == 0) % if its viable and the model predicts it won't grow
                % force it to grow
                disp('mismatch type NGG; editing network')
                lb11(find(weights)) = OPTIMAL_THRESH(ci)*f(ci); % it should atleast reach this value
                % run optimization
                [v00_1,f00_1,status1_1] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
                while (status1_1 ~= 5)
                    disp('no optimal solution found. relaxing expected growth rate by 1%')
                    % for tuning optimal_thresh.. if the expected growth is higher than the max. growth rate
                    % no optimal solution would be obtained. relax the predicted growth rate in that case
                    lb11(find(weights)) = lb11(find(weights)) * 0.99;
                    [v00_1,f00_1,status1_1] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
                    
                end
                disp('optimal thresh value used:')
                disp(lb11(find(weights))/f(ci)) % display the ratio that was finally used.
                
                % now check which flux changed among the targets
                temprxnpos1 = temprxnpos(constrainRxn); % these are the reactions for which i adjusted the vms
                if ~isempty(temprxnpos1)
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % debugging
                    size(temprxnpos1)
                    size(v00_1(temprxnpos1))
                    size(v00(ci,temprxnpos1)')
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        prob_flag = false;
                        clear spos sxpos ix_sx*

                    switch lower(metric_type)
                        case 'ratio threshold'
                            disp(metric_type)
			% ratio doesnt work coz of too many nans
%                             diff_ratio = abs ((v00_1(temprxnpos1) - v00(ci,temprxnpos1)')./v00(ci,temprxnpos1)');

                           diff_vec = v00_1(temprxnpos1) - v00(ci,temprxnpos1)';
% 		  use a z score based strategy for grouping.. ? 
			diff_ratio = abs(zscore(diff_vec));

			size(diff_ratio)

                            RATIO_THRESHOLD = 2 
                            spos = find(diff_ratio > RATIO_THRESHOLD);
                            % choose rxns at random from this list
                            spos = spos(randperm(length(spos)));
                                spos1 = find(diff_ratio <= RATIO_THRESHOLD);
                            spos1 = spos1(randperm(length(spos1)));
                            spos2 = [spos;spos1];
                            sxpos = temprxnpos1(spos2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basically spllitting the rxns into two groups; one that changes more and the other less; picking the more
% changing reactions at random first and then from the second section
                            % debugging
				disp(diff_ratio)
				disp(sxpos)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        case 'probability'
                            ix_sx1 = find(ismember(regulator,tfnames(tfstate)));
                            [sx, spos] = sort(probtfgene); 
                            sxpos = spos(ismember(spos,ix_sx1));
                            prob_flag = true;
                        otherwise
				disp('default mode for sorting reactions')                            
                            diff_vec = v00_1(temprxnpos1) - v00(ci,temprxnpos1)';
                            % sort by magnitude and find the rxns tat change most
                            [sxx, spos] = sort(abs(diff_vec),'descend'); 
                            % now remove this interaction
                            % and find the gene for the corresponding rxn
                            
                            sxpos = temprxnpos1(spos);
                    end
                    for remov = 1:length(sxpos)
                        
                        disp(['editing network. iteration no:  ',num2str(remov)]);
                        if prob_flag
                            ix_sx = sxpos(1:remov);
                            probtfgene(sxpos(1:remov)) = 1;
                            remove_interactions(sxpos(1:remov)) = 1;
                        else
                            gpos = genelist(ismember(rxnpos,sxpos(1:remov))); % these are the genes that are involved
                            % next i need to find the corresponding interactions with
                            % the tf
                            ix_sx1 = ismember(regulator,tfnames(tfstate));
                            ix_sx2 = ismember(regulated,model.genes(gpos));
                            ix_sx = (ix_sx1&ix_sx2);
                            probtfgene(ix_sx) = 1;
                            remove_interactions(ix_sx) = 1;
                        end
                        % now run PROM for the new network
                        reg1 = regulator(~remove_interactions);
                        reg2 = regulated(~remove_interactions);
                        ix_sx1 = ismember(reg1,tfnames(tfstate));
                        
                        if (remov < length(sxpos))&& (~isempty(reg1(ix_sx1)))
                            
                            [f_chk,f_ko_chk,v_chk,v_ko_chk,status1_chk,lostxns_chk] =  promv2(model,expression,expressionid,reg1(ix_sx1),reg2(ix_sx1),[],[],[],v11,v12,KAPPA,DATATHRESHVAL,[],sizeflag);
                            
                        else % all interactions are removed.. trivial solution.. doesnt happen usually
                            f_chk = f_wt;
                            remove_interactions1(ix_sx) = 1;
                            disp('trivial solution for TF no')
                            disp(ci)
                            break;
                        end
                        
                        % now check if it matches.. else loop around and do it again
                        
                        
                        if (f_chk >= OPTIMAL_THRESH(ci)*f_wt) % it worked
                            remove_interactions1(ix_sx) = 1;
                            disp('modified growth rate')
                            f(ci) = f_chk;
                            disp(f(ci))
                            break;
                        end  % if that didnt work i'll run the for loop adding the next interaction
                        
                    end
                    
                    
                end
            end
        else
            disp('predicted phenotype matches measured phenotype')
        end
        
    else
        disp('phenotype data not available for comparison')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    coun  = 1;        lbh = lbg; ubh  = ubg;
    
    while ((status1(ci) ~= 5) || (v00(ci,find(weights)) < 0))
        lbh(lbh ~= lbf) = lbh(lbh ~= lbf) - 1E-3;
        ubh(ubh ~= ubf) = ubh(ubh ~= ubf) + 1E-3;
        dxdt0 = [zeros(size(S,1),1);lbh;ubh];
        [v00(ci,:),f00(ci),status1(ci)] = glpk(-weights11,A,dxdt0,lb11,ub11,ctype1,[],[],param);
        coun = coun + 1
        if (coun > 2),
            % if its impossible to estimate, then
            % check the unweighted growth rate and the one with max weight prom -
            % if very less difference use it - any way warn the user about
            % the problem at the iteration number - ci;
            [v3,f3,status3] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
            [v30,f30,status30] = glpk(-weights00,A,dxdt0,lb11,ub11,ctype1);
            %if abs((f3- v30(find(weights)))/abs(f3)) < 0.1
            f00(ci) = -f3;
            v00(ci,find(weights)) = -f3;
            %else
            disp(' problem running PROM in'); disp(ci);break;  % if that doesnt work,  display a warning
            %end
            
            disp('check network for TF no '); disp(ci); break;
            
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lbg_st(ci,:) = lbg;
    ubg_st(ci,:) = ubg;
    
    lb_st(ci,:) = lbff;
    ub_st(ci,:) = ubff;
    
    [v2(ci,:),f1(ci),status] = glpk(-weights,S,dxdt,lbg,ubg,ctype);
    f_ko(ci) = -f1(ci);
    
    ktime = toc;
    waitpar = [num2str(ceil(ci/length(tfnames)*100)),'% complete. Time taken:',num2str(ceil(ktime)),' secs'];
    waitbar(ci/length(tfnames),hw,waitpar);
    
    
    ff00(scou,ci) = v00(ci,find(weights));
    %disp(ff00(scou,ci))
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if datathreshflag
        if all(lost_xn(k))   % if none of the probabilities of a gene can be estimated, then ->
            disp('Error: interactions cant be estimated; change binarization threshold')
            v00(ci,:) = NaN;
            f1(ci) = NaN;
            v2(ci,:) = NaN;
            f00(ci) = NaN;
            %break;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear tempgenepos tempgeneprobs temprxnpos k
    %disp(bnumstobekoed(ci))
end

% f_ko = -f1';
v_ko = v2;
%
% f00_ko(:,scou) = v00(:,find(weights));
v00_ko = v00;
%
% f = f00_ko;
v = v00_ko;


lostxns(:,scou) = lost_xn;

%disp([f,f_ko])

end

