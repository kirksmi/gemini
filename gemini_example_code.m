%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example code for running gemini on all yeast TFs and
% glucose minimal media
% simple version
load yeast_gemini_data_sc1 regulator targets expression expressionid model v11 v12 phenotype  
[f,inital_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype,'lethal',[],v11,v12,[],0.05,[],0.33); % default params - 0.33 for binarizing expression; 0.05 for calling lethal/non lethal phenotypes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example code for running gemini on a small subset of yeast TFs and
% glucose minimal media
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load yeast_gemini_data_sc1 regulator targets expression expressionid model v11 v12 phenotype  
te = unique(regulator);
% running gemini for a small subset - may be on a cluster
li = 0:ceil(length(te)/40):length(te);
li(end+1) = length(te);
split = 1; 
indi = li(split)+ 1:li(split+1)
st = clock;
sd = [date,'_',num2str(st(4))];
subsets = te(indi);
targets = targets(ismember(regulator,subsets));
    regulator = regulator(ismember(regulator,subsets));
phenotype = phenotype(indi);
%kappavec = [0.01,0.1,1,10,100,1000,10000]; 
%optthreshvec = [0.03;0.04;0.05;0.06;0.07;0.08];
%datathreshvals = [0.1,0.2,0.25,0.33,0.4,0.5,0.75,0.9];

% for doing flux variability
if (length(regulator) < 1000) % for small networks its faster to compute min/max for the subset thats regulated than the entire metab model
    sizeflag = 1;
else
    sizeflag = 0;
end


   
   
[f,inital_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype,'lethal',[],v11,v12,sizeflag,0.05,[],0.33); % default params - 0.33 for binarizing expression; 0.05 for calling lethal/non lethal phenotypes

   tempname = ['save gemini_results_',sd];


eval(tempname)



