# Integrating genetic regulatory networks and genome-scale metabolic models using GEMINI and PROM
We have developed two genome-scale metabolic network algorithms that integrate the transcriptional regulatory network into genome-scale metabolic models. 

## PROM (Probabilistic Regulation of Metabolism) 
**PROM** enables the quantitative integration of regulatory and metabolic networks to build genome-scale integrated metabolic–regulatory models. 

## Installation
PROM was implemented in MATLAB (recommended version: 2018+), and requires the GLPK Solver (recommended version: 2018+). The file `promv2.m` in this repository is used to run PROM.

## Usage
The function `promv2` has the following syntax:
`function [f,f_ko,v,v_ko,status1,lostXns,probTFgene] =  promv2(model,expression,expressionID,regulator,targets,litEvidence,prob_prior,subsets,minFlux,maxFlux,KAPPA,DATATHRESHVAL,probTFgene,sizeFlag)`
#### INPUTS:
* `model`: The metabolic model obtained from COBRA toolbox through the readcbmodel command. 

* `expression`: Gene expression data. 
  * Rows are genes, columns are conditions.
  * There is no need to normalize or impute.

* `expressionID`: An array of identifiers for each row/gene.

* `regulator`, `targets`: Regulatory network.
  * Format: cell array of regulators and matching target genes.
  * Example:   
           `Regulator = {'RegA'; 'RegB' ; 'RegC'};`  
           `Targets = {'GeneA';'GeneB';'GeneC'};`
  * Note that the names or identifiers used in the regulatory data should match the names/IDs given for gene expression data.

 * `minFlux`, `maxFlux`: Obtained from either `fastfva` or `fluxvariability` command in COBRA.
   * `[minFlux, maxFlux] = fastFVA(model);`

* `sizeFlag`: Tells PROM if the regulatory network is large. 
  * 0 for large networks
  * 1 for small networks ( less than 1000 interactions)

OPTIONAL PARAMETERS

* `litEvidence`: High confidence interactions (not necessarily based on literature) should be flagged as 1 in `litevidence` array.
  * The remaining interactions should be set to 0.
  * Should have same length as the regulator/target arrays.
  
* `probPrior`: Should be set to values between 0 and 1 for those interactions with litEvidence (other values in the array will be ignored)
  * Should have same length as the regulator/target arrays.

* `KAPPA`: determines strength of regulatory constraints.
  * Default value is 1, which should work for most systems.
  
* `subsets`: Subsets of TFs for which PROM should be run.
  * Default: run for all TFs.

#### OUTPUTS
* `f`: Growth rate.

* `f_KO`: Growth rate after regulator knockout.

* `v`: Flux response.

* `v_KO`: Flux response after regulator knockout.

* `status1`: The glpk solver status.
  * The status should be 5 for glpk.

* `lostXns`: Gives the interactions that could not be quantified based on the
threshold set for binarization.
  * The program will output a warning if the threshold is poorly chosen. The default value (0.2 - 0.4) should work for most cases.

* `probTFgene`: Gives the probabilities estimated for each interaction.


## GEMINI (Gene Expression and Metabolism Integrated for Network Inference)
**GEMINI** directly connects regulatory interactions to observable phenotypes and allows rapid assessment of inferred regulatory interactions using a metabolic network.

## Contributions
Contributions are welcome! Please read the contributions guide to get started. Also feel free to submit bugs, feature requests, and pull requests.

Additionally, you can support development for PROM by citing the original publications.

## Publications
1. Chandrasekaran S and N.D. Price, "Probabilistic integrative modeling of genome-scale metabolic and regulatory networks in Escherichia coli and Mycobacterium tuberculosis," PNAS, 2010. 
2. Simeonidis E, Chandrasekaran S, Price ND. “A guide to integrating transcriptional regulatory and metabolic networks using PROM (Probabilistic Regulation of Metabolism)”, Methods in Molecular Biology: Systems Metabolic Engineering.
3. Chandrasekaran S and N.D. Price, “Metabolic Constraint-based Refinement of Transcriptional Regulatory Networks”, PLOS Computational Biology, 2013.
