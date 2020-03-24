GEMINI 
    
GEMINI (Gene Expression and Metabolism Integrated for Network Inference) produces a regulatory network that is simultaneously consistent with observed gene knockout phenotypes, gene expression data, and the corresponding metabolic network.

~~~~~~~~~~~
 >> [f,initial_network,final_network] =  GEMINI(model,expression,expressionid,regulator,targets,phenotype,type,subsets,v11,v12,sizeflag,OPTIMAL_THRESH,metric_type)

 GEMINI call with minimal number of inputs
>> [f,initial_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype);

~~~~~~~~~~~~

 INPUTS

 Model is the metabolic model for the organism (obtained from COBRA toolbox through readcbmodel command) . The model should be set to a specific growth condition under study ( like glucose minimal media)

 Gene expression data - rows - genes,columns - conditions; (preferably normalized and imputed)

 Expressionid - an array of identifiers for each row/gene should be included

draft regulatory network - format - cell array of regulators and matching target genes.
For example,   
~~~~~~~~~~~~~~~~~~~
>> Regulator = {'RegA'; 'RegB' ; 'RegC'};  Targets ={'GeneA';'GeneB';'GeneC'};
~~~~~~~~~~~~~~~~~~~~
note that the names or identifiers used in the regulatory data should match the names/ids given for gene expression data

 phenotype - logical vector (true/false) - the growth phenotype of each transcription factor knockout under a specific condition specified by the metabolic model. it should be the same length as the number of TFs in the model. 


OPTIONAL

type - a string describing the knockout phenotype data - should be either 'lethal' (default) or 'suboptimal' ; this is required to set the OPTIMAL_THRESH value

OPTIMAL_THRESH is the threshold for determining lethal/non-lethal phenotypes; by default, a knockout that grows less than 5% of the wild type growth rate is considered lethal and less than 95% of the wildtype growth rate is considered suboptimal; ( OPTIMAL_THRESH default value = 0.05 for 'lethal' type, corresponds to 5% and 0.95 for 'suboptimal')

subsets : subsets of tfs for which GEMINI should be run ; default - for all tfs

v11,v12 are the minimum and maximum possible flux through each reaction in the model for the given condition. this obtained through flux variability analysis (either from fastfva or fluxvariability command in COBRA toolbox)

~~~~~~~~~~~~~~~ 
>> [v11,v12] = fluxvariability(model); 
~~~~~~~~~~~~~~~~

 OUTPUT 

the algorithm gives the refined network and the growth rate (f) after knock out of all tfs in the refined network.  note that growth rate(f) is only semi-quantitative, unless trained on suboptimal growth data

 EXAMPLE

~~~~~~~~~~~~~~~~~~~ 
>> load yeast_gemini_data_sc1 regulator targets expression expressionid model v11 v12 phenotype  

>> [f,initial_network,final_network] = GEMINI(model,expression,expressionid,regulator,targets,phenotype,'lethal',{'YAL051W'},v11,v12);
~~~~~~~~~~~~~~~~~~~~~~~~~