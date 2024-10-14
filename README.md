# PAN-GO_CAFA_evaluation
Contains code and data used to perform CAFA evaluation of PAN-GO data

Annotations from six different automated annotation methods were compared with those from Pan-GO.
* DeepGOCNN
* DeepGOPlus
* DeepGOZero
* Diamond
* PANNZER
* TALE


## A. Input files:

1. `PanGO_curated_genes` - This file contains a list of human genes in PanGO families that have at least one experimental annotation to any gene in the family. This is considered a list of curatable or trainable human genes. The file is also used to map the UniProt ID to the PANTHER gene IDs.
2. `training_exp_human` - All experimental annotations to human genes from Gene Ontology [October 2019 release](https://release.geneontology.org/2019-10-07/index.html). The file includes directly annotated terms, as well as all the inferred parent terms based on GO ontology file. When parent terms are inferred, only 'is_a' and 'part_of' relationships were used.
These annotations, together with experimental annotations to all genomes, were used for training.
The following terms (root terms and some binding terms) were excluded from the file: 
    * GO:0005515	protein binding
    * GO:0005488	binding
    * GO:0003674	molecular function
    * GO:0008150	biological process
    * GO:0005575	cellular component
3. `new_annotations` - These are new experimental annotations from the [March 2022 GO release](https://release.geneontology.org/2022-03-22/index.html). Annotations to the parent terms were inferred based on GO hierarchy.
Annotations, including those to the parent GO terms, were excluded if they were already in the `training_exp_human` file.
The following terms (root terms and some binding terms) were excluded from the file: 
    * GO:0005515	protein binding
    * GO:0005488	binding
    * GO:0003674	molecular function
    * GO:0008150	biological process
    * GO:0005575	cellular component

These annotations are also known as "ground truth" in CAFA.

4. `go_all_parents` - A lookup file of all parent terms for each GO term. The file was generated based on the [go.obo file](https://release.geneontology.org/2022-03-22/ontology/go.obo) released in March 2022.

5. Prediction files (`prediction_[method]`) - For the six methods to be compared, the predictions were made according to the instructions from each method. Experimental annotations from the GO October 2019 release were used for training.
Annotations, including those to the parent GO terms, were excluded if they were already in the `training_exp_human` file.
The following terms (root terms and some binding terms) were excluded from the prediction files: 
    * GO:0005515	protein binding
    * GO:0005488	binding
    * GO:0003674	molecular function
    * GO:0008150	biological process
    * GO:0005575	cellular component

## B. Scripts

1. `PanGO_CAFA_F_score.pl` - A script to generate precision, recall and f_score for each prediction method.
2. `PanGO_comparison.pl` - A script to compare annotations between Pan-GO and each automated prediction method.

## C. Usage
```
perl PanGO_CAFA_F_score.pl -i prediction_[method] -t training_exp_human -n new_annotations -g PanGO_curated_genes -p go_all_parents -o [method] > output_[method]

perl PanGO_comparison.pl -i prediction_PanGO -g PanGO_curated_genes -p go_all_parents -a prediction_[method] > [output_file]
```
