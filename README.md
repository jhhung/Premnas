# Premnas
#### Premnas is a Framework to Estimate Intraclonal Heterogeneity and Subpopulation Changes from Perturbational Bulk Gene Expression Profiles.

## Information

Workflow of Premnas: 

![](https://i.imgur.com/sLydog1.png)


Three main steps are included in Premnas:
1. Learning ad hoc subpopulation characteristic
2. Performing digital cytometry
3. Analyzing subpopulation change

To execute Premnas, we provide a script for Premnas users. The script would automatically downlaod two required dockers and run through the whole Premans' workflow.

## Requirement
* [Docker](https://www.docker.com/) version higher than 19.03.8
* [Python 3.7 ](https://www.python.org/downloads/) or higher as well as [numpy](https://numpy.org/) and [pandas](https://pandas.pydata.org/) packages
* Registration on [CIBERSORTx](https://cibersortx.stanford.edu/index.php) and get token for CIBERSORTx docker container on https://cibersortx.stanford.edu/getoken.php .



## Execution

```sh
> git clone https://github.com/jhhung/Premnas.git
> cd Premnas
    
#run with the shell script:
> sh Premnas.sh \
    -I {absolute_input_dir_path} \
    -O {absolute_output_dir_path} \
    -D {single_cell_GEPs_data} \
    -S {single_cell_source_clone} \
    -U {CIBERSORTx_username} \
    -T {CIBERSORTx_token} \
    -C {CMap_data} \
    -M {CMap_metadata} \
```

Parameter explaination:
* Required

    * -I [absolute_input_dir_path]:
        Absolute folder path to your all input files, such as single cell GEPs and CMap data.
        
    * -O [absolute_output_dir_path]: 
        Absolute folder path to your all output files.
        
    * -D [Single_cell_GEPs_data]: 
        Single cell GEPs data used to learn the subpopulation characteristic.
        
    * -S [single_cell_source_clone]: 
        Source clone labels for all cells. In Premnas, we perform batch corretion for single cells by their source clone.
        
    * -U [CIBERSORTx_username]:
        Registered user name of [CIBERSORTx](https://cibersortx.stanford.edu/index.php). 
        
    * -T [CIBERSORTx_token]: 
        According to CIBERSORTx team policy, it is necessary to run CIBERSORTx docker with an authentication token. The token can be apply on https://cibersortx.stanford.edu/getoken.php .
        
    * -C [CMap_data]: 
        LINCS L1000 CMap level 3 gene expression profiles which can be downloaded from [the GEO website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) or other bulk gene expression profiles. 
        
    * -M [CMap_metadata]: 
        Information of each perturbation experiment in the LINCS L1000 CMap database (e.g., 'GSE70138_Broad_LINCS_inst_info_2017-03-06.txt.gz' on [the GEO website](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138))

* Optional.

    * -m [Learning_Mode] (default: false): 
        If set true, the script would only learn ad hoc subpopulation characteristic but not perform digital cytometry. The learning mode is designed for users who want to do digital cytometry by other tools. Note that you could ommit parameters -U, -T, -C, -M in learning mode.
        
    * -s [Susceptibility_threshold] (default: 0.9): Threshold used to evaluate the inhibitory effects of each perturbagen-concentration-time pair (PCT pair) for each subpopulation. The value should be set between 0 and 1.
 
    * -c [Consistency_threshold] (default: 0.8): Threshold used to check whether the perturbagen of the PCT pair with higher doses would also perform well. Each PCT pair should past this consistency test before adding to cocktail therapy.

Premnas would generate several files under the ```absolute_output_dir_path```, which user has set previously:
```c
output_dir   
    |--- Subpopulation-charactistic.txt 
    |       // Main output, would be used for digital cytometry (CIBERSORTx)
    |--- ACTIONet-model
    |       // Model learn by ACTIONet
    |--- assigned-subpopulation.txt
    |       // Subpopulation labels of each cells assigned by ACTIONet
    |--- pruned-assigned-subpopulation.txt
    |       // Subpopulation labels of each cells after pruning some cells by considering archetypal explicit function
    |--- subpopulation-visualization.pdf
    |        // UMAP visualization of single cells.
    |--- CIBERSORTx_Adjusted.txt
    |        //Deconvolution result generate by CIBERSORTx
    |--- Treatment-selection-output.csv
            //Eventually selected cocktail therapy
            
```

## Reference
1. LINCS L1000 CMap: Subramanian, Aravind et al. A Next Generation Connectivity Map: L1000 Platform and the First 1,000,000 Profiles. Cell vol. 171,6 : 1437-1452.e17. doi:10.1016/j.cell.2017.10.049 (2017)
2. ACTIONet: Mohammadi, S., Davila-Velderrain, J. & Kellis, M. A multiresolution framework to characterize single-cell state landscapes. Nature Communications 11, 5399, doi:10.1038/s41467-020-18416-6 (2020).
3. CIBERSORTx: Newman, A. M. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37, 773-782, doi:10.1038/s41587-019-0114-2 (2019).
