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
* [cmapPy](https://pypi.org/project/cmapPy/) should be installed if expect CMap database automatically downlaod. 



## Execution

User could decide to execute the whole workflow of Premnas or just run the specific step by selecting the "mode" argument, which are respectively descripted in more detail below. 

Moreover, to make the data input meet our file format, user could also download CMap database with specific cell line by our script. The method of the automatically download mode is also descripted below.

### Run the whole Premnas framework

```
> git clone https://github.com/jhhung/Premnas.git
> cd Premnas
    
> python3 Premnas.py All \
    -input_dir {absolute_input_dir_path} \
    -output_dir {absolute_output_dir_path} \
    -single_cell {single_cell_GEPs_data} \
    -sc_source {single_cell_source_clone} \
    -user_name {CIBERSORTx_username} \
    -token {CIBERSORTx_token} \
    -mixture {CMap_data} \
    -sub_characteristic {subpopulation_characteristic_file}
```

The description of each argument could be founded below the README file or by `python3 Premnas.py -help`. 

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

### Run the specific step of Premans

#### 1. Learning ad hoc subpopulation characteristic:

```
> python3 Premnas.py Learn_char \
    -input_dir {absolute_input_dir_path} \
    -output_dir {absolute_output_dir_path} \
    -single_cell {single_cell_GEPs_data} \
    -sc_source {single_cell_source_clone} 
```
Parameter explanation:
* -input_dir [absolute_input_dir_path]:
        Absolute folder path to your all input files, such as single cell GEPs and CMap data.
        
* -output_dir [absolute_output_dir_path]: 
        Absolute folder path to your all output files.   

* -single_cell [Single_cell_GEPs_data]: 
        Single cell GEPs data used to learn the subpopulation characteristic.

* -sc_source [single_cell_source_clone]: 
        Source clone labels for all cells. In Premnas, we perform batch corretion for single cells by their source clone.

#### 2. Performing digital cytometry

```
> python3 Premnas.py Dig_cytometry \
    -input_dir {absolute_input_dir_path} \
    -output_dir {absolute_output_dir_path} \
    -user_name {CIBERSORTx_username} \
    -token {CIBERSORTx_token} \
    -mixture {CMap_data} \
    -sub_characteristic {subpopulation_characteristic_file}
```

Parameter explanation:
* -input_dir [absolute_input_dir_path]:
        Absolute folder path to your all input files, such as single cell GEPs and CMap data.
        
* -output_dir [absolute_output_dir_path]: 
        Absolute folder path to your all output files.   

* -user_name [CIBERSORTx_username]:
        Registered user name of [CIBERSORTx](https://cibersortx.stanford.edu/index.php). 

* -token [CIBERSORTx_token]: 
        According to CIBERSORTx team policy, it is necessary to run CIBERSORTx docker with an authentication token. The token can be apply on https://cibersortx.stanford.edu/getoken.php .

* -mixture [CMap_metadata]: 
        CMap metadata used to map probe IDs to gene names.
        
* -sub_characteristic [Subpopulation characteristic]: 
        Subpopulation characteristic file, which could generate by the previous step of Premnas.

#### 3. Analyzing subpopulation change

```
> python3 Premnas.py Analyze_sub \
    -input_dir {absolute_input_dir_path} \
    -output_dir {absolute_output_dir_path} \
    -metadata {CMap_metadata}
```
Parameter explanation:
* -input_dir [absolute_input_dir_path]:
        Absolute folder path to your all input files, such as single cell GEPs and CMap data.
        
* -output_dir [absolute_output_dir_path]: 
        Absolute folder path to your all output files.
        
* -metadata [CMap_metadata]: 
        CMap metadata used to map probe IDs to gene names.

### Download the CMap database with specific cell line

```
> python3 Premnas.py Download_CMap_data \
    -celltype {Cell_line_argument}
```
The `-celltype` could be set as:
* A375
* A549
* HCC515
* HEPG2
* MCF7
* PC3
* VCAP
* HT29

The python script would download LINC L1000 CMap data with specific cell line including its metadata (which would be used in "Analyzing subpopulation change" step in Premnas) from GEO website.   

        

## Reference
1. CMap: Lamb, J. et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science 313, 1929-1935, doi:10.1126/science.1132939 (2006).

2. ACTIONet: Mohammadi, S., Davila-Velderrain, J. & Kellis, M. A multiresolution framework to characterize single-cell state landscapes. Nature Communications 11, 5399, doi:10.1038/s41467-020-18416-6 (2020).
3. CIBERSORTx: Newman, A. M. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37, 773-782, doi:10.1038/s41587-019-0114-2 (2019).