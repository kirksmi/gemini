# Integrating genetic regulatory networks and genome-scale metabolic models using GEMINI and PROM
We have developed two genome-scale metabolic network algorithms that integrate the transcriptional regulatory network into genome-scale metabolic models. 

**PROM (Probabilistic Regulation of Metabolism)** enables the quantitative integration of regulatory and metabolic networks to build genome-scale integrated metabolic–regulatory models. 

## Installation
PROM was implemented in MATLAB (recommended version: 2018+), and requires the GLPK Solver (recommended version: 2018+). The file `promv2.m` in this repository is used to run PROM.

## Usage
The function `promv2` has the following syntax:
`function [f,f_ko,v,v_ko,status1,lostxns,probtfgene] =  promv2(model,expression,expressionid,regulator,targets,litevidence,prob_prior,subsets,v11,v12,KAPPA,DATATHRESHVAL,probtfgene,sizeflag)`
### INPUTS:
* `model`: - The metabolic model obtained from COBRA toolbox through the `readcbmodel` command.

* `expression`: Gene expression data.
  * Rows are genes, columns are conditions.
  * There is no need to normalize or impute.

* `expressionid`: - An array of identifiers for each row/gene.

* `regulator`, `targets`: Regulatory network
  * Format: cell array of regulators and matching target genes.
  * Example:   
           `Regulator = {'RegA'; 'RegB' ; 'RegC'};`  
           `Targets = {'GeneA';'GeneB';'GeneC'};`
  * Note that the names or identifiers used in the regulatory data should match the names/IDs given for gene expression data.

 * `v11`, `v12`: Obtained from either `fastfva` or `fluxvariability` command in COBRA.
   * `[v11, v12] = fastFVA(model);`

* `sizeflag`: Tells PROM if the regulatory network is large
  * 0 for large networks
  * 1 for small networks ( less than 1000 interactions)

OPTIONAL PARAMETERS

* `litevidence`: High confidence interactions (not necessarily based on literature) should be flagged as 1 in `litevidence` array.
  * The remaining interactions should be set to 0.
  * Should have same length as the regulator/target arrays.
  
* `prob_prior`: Should be set to values between 0 and 1 for those interactions with litevidence (other values in the array will be ignored)
  * Should have same length as the regulator/target arrays.

* `KAPPA`: determines strength of regulatory constraints.
  * Default value is 1, which should work for most systems.
  
* `subsets`: Subsets of TFs for which PROM should be run.
  * Default: run for all TFs.

### OUTPUTS
* The algorithm gives the **growth rate** `f` and **flux response** `v` after knock out of all regulators in the regulatory model.
* `status`: The glpk solver status.
  * The status should be 5 for glpk.

* `lostxns`: Gives the interactions that could not be quantified based on the
threshold set for binarization.
  * The program will output a warning if the threshold is poorly chosen. The default value (0.2 - 0.4) should work for most cases.

* `probtfgene`: Gives the probabilities estimated for each interaction.



**GEMINI (Gene Expression and Metabolism Integrated for Network Inference)** directly connects regulatory interactions to observable phenotypes and allows rapid assessment of inferred regulatory interactions using a metabolic network.

### Publications
1. Chandrasekaran S and N.D. Price, "Probabilistic integrative modeling of genome-scale metabolic and regulatory networks in Escherichia coli and Mycobacterium tuberculosis," PNAS, 2010. 
2. Simeonidis E, Chandrasekaran S, Price ND. “A guide to integrating transcriptional regulatory and metabolic networks using PROM (Probabilistic Regulation of Metabolism)”, Methods in Molecular Biology: Systems Metabolic Engineering.
3. Chandrasekaran S and N.D. Price, “Metabolic Constraint-based Refinement of Transcriptional Regulatory Networks”, PLOS Computational Biology, 2013.
