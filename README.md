# Premnas
#### Premnas is a Framework to Estimate Intraclonal Heterogeneity and Subpopulation Changes from Perturbational Bulk Gene Expression Profiles.

## Information

Workflow of Premnas: 
![](https://i.imgur.com/D0Nacge.png)

Three main steps are included in Premnas:
1. Learing ad hoc subpopulation characteristic
2. Performing digital cytometry
3. Analyzing subpopulation change

To execute the whole Premnas, we provide the dockers and source code for Premnas users. The way to run above three steps are introduced below respectively. 


## Execution

### Learing ad hoc subpopulation characteristic
    
* Run with docker
    ```cpp=
    docker pull JHHLAB/Premnas-learn-sub
    docker run -v {input_dir_path}:/input_dir \-v {output_dir_path}:/output_dir JHHLAB/Premnas-learn-sub
    ```
    Make sure there are two files under your input directory path: 
    * ```Single-cell-data.csv``` : M X N single cell data with M represents the number of genes and N represents the number of cells. 
    * ```Metadata.txt``` : N X 2 matrix, which each row includes a cell name and the correesponding source of clone.

    The docker would automatically reads above files as input by their file name.
    
* Run with script
    ```python=
    git clone https://github.com/jhhung/Premnas.git
    cd Premnas
    chmod 777 Premnas.sh {Single-cell-data.csv} {Metadata.txt}
    ./Premnas.sh
    ```
    
The outputs of the docker would be like:
```c=
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
            // UMAP visualization of single cells.
    
```

### Digital cytometry

Premnas perform digital cytometry through CIBERSORTx. The CIBERSORTx team has provided docker, which could download in https://cibersortx.stanford.edu/download.php, to allow user to execute CIBERSORTx.

In Premnas, we run the CIBERSORTx docker by:

```c=
docker pull cibersortx/fractions
docker run \
    -v {input_dir_path}:/src/data \
    -v {output_dir_path}:/src/outdir \  
     cibersortx/fractions \
    --username email_address_registered_on_CIBERSORTx_website \
    --token token_obtained_from_CIBERSORTx_website \
    --single_cell TRUE \
    --refsample {Single-cell-profile.txt} \ 
    --mixture {Bulk.txt} --fraction 0 \
    --rmbatchSmode TRUE \
    --perm 500
```

Note that the ```Single-cell-profile.txt``` should be the ```Subpopulation-charactistic.txt ``` which has generated in previous step, and the ```bulk.txt``` should be the CMap database you want to analysis.

The CIBERSORTx docker would generat a series of file including ```CIBERSORTx_Adjusted.txt```, which represented the deconvolution results of CMap data and would further be used in next step.

### Analyzing subpopulation change

```python=
git clone https://github.com/jhhung/Premnas.git
Treatment-selection.py \
    {CMap-metadata.txt} \
    CIBERSORTx_Adjusted.txt 
```
* ```CMap-metadata.txt``` should be the metadata provided in CMap database. It records the perturbagen dose and time for each sample in CMap. Note that you should input the right metadata, which is corresponding to the CMap database you analysis previous digital cytometry step.

The eventually cocktail therapy would be store in file ``` Treatment-selection-output.csv```.

## Reference
1. CMap: Lamb, J. et al. The Connectivity Map: using gene-expression signatures to connect small molecules, genes, and disease. Science 313, 1929-1935, doi:10.1126/science.1132939 (2006).

2. ACTIONet: Mohammadi, S., Davila-Velderrain, J. & Kellis, M. A multiresolution framework to characterize single-cell state landscapes. Nature Communications 11, 5399, doi:10.1038/s41467-020-18416-6 (2020).
3. CIBERSORTx: Newman, A. M. et al. Determining cell type abundance and expression from bulk tissues with digital cytometry. Nature Biotechnology 37, 773-782, doi:10.1038/s41587-019-0114-2 (2019).



