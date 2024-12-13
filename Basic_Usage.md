# **Basic Usage**

<p align='justify'>
This tutorial demonstrates how to use MAYA with a dataset annotated with SMILES notation and associated activity or property data to automate the analysis of SMARts patterns and generate a chemical multiverse. The purpose of this tutorial is to serve as a learning tool for using MAYA. To enhance accessibility and efficiency for the user, we simplify the code analysis and focus on practical implementation.

### Set up compounds dataset  
<p align='justify'>
The current version of **MAYA** supports several file types, including .CSV, .XLSX, .TSV, .XLSX, .JSON and .XML.  
The dataset must be annotated with SMILES notation and include at least 1 activities or properties. The structure of the input is illustrated in the dataset examples. We need to define certain variables to ensure that MAYA uses the correct columns for calculations.

1. **dataset** (CSV, XLSX, TSV, JSON, XML): DataFrame of compounds provided by the user.

2. **ID** (str): Column name containing the ID of each compound. 

3. **smiles\_column\_name** (str): Column name containing SMILES notation.

4. **target\_activities** (list): List with the names of columns containing target activity values.

5. **vPCA** (bool): User-defined variable for applying Principal Component Analysis., defined True as default

6. **t\_SNE** (bool): User-defined variable for applying t-SNE, defined True as default

7. **MACCS** (bool): User-defined variable for applying MACCS keys calculation, defined True as default

8. **ECFP** (bool): User-defined variable for applying Extended Connectivity Fingerprint calculation., defined True as default

9. **radius** (int): User-defined variable for ECFP radius 2 or 3., defined 3 as default

10. **fpsize** (int): Variable for define the size of ECFP, the number of bits, defines 2048 as default

11. **druglikeness\_descriptors** (bool): User-defined variable for applying the 6 drug-likeness descriptors, defined True as default

12. **signaturizer\_code** (list): User-defined variable for applying signaturizer with 6 different codes.

13. **palette** (str): User-defined variable for defining a continuous color palette., ‘RdBu\_r’ as default

14. **perplexity** (int): User-defined variable for adjusting the proximity of points in t-SNE, depending on dataset size, defined 33 as default

15. **size\_point** (float): User-defined variable for adjusting point size in the plot, defined 13 as default

16. **point\_shape** (str): User-defined variable for adjusting point shape in the plot, ‘circle’ as default

## Relevance of some variables:
<p align='justify'>
   
##### **radius**  
<p align='justify'>

Defined as the radius of the circle centered on each atom within a molecule. Typically this radius corresponds to 2 or 3 bonds, allowing the mapping of the immediate chemical environment of each atom. This variable is essential for capturing local chemical interactions that contribute to molecular activity.

#### **size_point_representation**
<p align='justify'>
The default point size is 13. However, it is possible to adjust the size to represent a specific characteristic of the database. By default, the size will represent the standard deviation of activity values, indicating the variability of activity values across the multiple targets in the databases. The user can modify the characteristic to be represented by defining the variable size_point_representation. The possible representations are as follows:

1. normal_desviation: Default configuration
2. size_point_HBA: Hydrogen Bond Acceptor
3. size_point_LogP: Octano - Water Coefficient
4. size_point_TPSA: Total Polar Surface Area
5. size_point_MW: Molecular Weight
6. size_point_HBD: Hydrogen Bond Donor
7. size_point_RB: Rotable Bonds

#### **signaturizer\_code**  
A list containing the names of the different bioactive descriptors[^1] used. These descriptors can be expanded or reduced depending on the analysis requirements. The included descriptors are:

a. Level A: Chemistry  
These descriptors are derived from the chemical structure

A1: 2D fingerprints  – A2: 3D fingerprints  – A3: Scaffolds  – A4: Structural keys  –  A5: Physicochemistry

b. Level B: Targets  
Descriptors related to the interactions of a compound with its biomolecular target, including data from biochemical assays and databases such as ChEMBL, PubChem or BindingDB.  
   
B1: Mechanisms of action  – B2: Metabolic genes – B3: Crystals – B4: Binding – B5: HTS bioassays

c. Level C: Biological Networks  
Descriptors designed to characterize how compounds modulate or affect interconnected biological systems, such as metabolic networks

C1: Small molecule roles – C2: Small molecule pathways – C3: Signaling pathways – C4: Biological processes

d. Level D: Cells  
Describes the impact of bioactive molecules on cellular processes

D1: Transcription – D2: Cancer lines – D3: Chemical genetics – D4: Morphology – D5: Cell bioassays

e. Level E: Clinic  
Describe the effects observed in human or animal models, including both therapeutic and side effects

E1: Therapeutic areas – E2: Indications – E3: Side effects – E4: Diseases & toxicology – E5: Drug-Drug interactions

#### **palette**  
This variable is used to modify the color palette applied in the graphs. It supports two main configurations\>

1. Predefined Matplotlib palettes: Users can select one of the color palettes already available in the Matplotlib library, such as ‘viridis’, ‘plasma’ or ‘coolwarm’.  
2. Custom list: Alternatively, a custom list of colors can be defined using RGB or HEX codes (e.g. \#FF5733 for orange).

This flexibility allows users to tailor visualizations to specific needs or aesthetic requirements.

#### **perplexity**  
The perplexity value is closely related to the resulting visualizations. This parameter represents an approximation of the number of nearest neighbors for each point, balancing the weight between local and global features in the analysis. It is recommended to set a perplexity value between 5 and 50 depending on the database size, ensuring that the value is always less than the total size of the database

## Scale
The visual representation of a property or activity is implemented to facilitate the identification of relevant patterns, differences, or trends. By default, the color scale in the graphs represents the sum of the pIC50 values if each compound's activities, defined as:
<p align='center'>
mpIC5o = -log(IC50(A)) + -log(IC50(A)) + ...
<p align='justify'>
However, the user can modify this configuration to visualizae the following properties:
   
1. normal_desviation: Standard deviation of each compounds's activities
2. HBD: Hydrogen Bond Donors
3. HBA: Hydrogen Bond Acceptors
4. MW:Molecular Weight
5. TPSA: Total Polar Surface Area
6. Logp: Octanol/Water Partition Coefficient
7. RB: Rotable bonds
8. mpIC50: Default configuration

Additionally, weighting and penalizations considerations can be implemented to enhance the analysis copabilities based on user requeriments. This allows the generation of a ranking or evaluation matric based on the activities or properties available in the database. The **evaluation_metric** and **metric_name** variables have been defined.

Defining a custom metric: Users can specify how to combine activities or properties. For example, for a series of three activities (Activity_1_Colum_Name, Activity_2_Colum_Name, Activity_3_Colum_Name), where one is undesired (Activity_3_Colum_Name), it is possible to penalize it, defining **evaluation_metric**. For exmaple:

```markdown
# To penalize an undesired activity
evaluation_metric = ('Activity_1_Colum_Name') + ('Activity_2_Colum_Name') - ('Activity_1_Colum_Name')
```

```markdown
# To assign different weights to activities/properties
evaluation_metric = 2 * ('Activity_1_Colum_Name') + ('Activity__Colum_Name') - ('Activity_3_Colum_Name')
```

```markdown
# To work with logaritmic values:
evaluation_metric = -log('Activity_1_Colum_Name') + 0.5(-log('Activity_2_Colum_Name')) - (-log('Activity_3_Colum_Name'))
evaluation_metric = -log('Activity_1_Colum_Name' * 1000) + (-log('Activity_2_Colum_Name' * 1000)) - (-log('Activity_3_Colum_Name' * 1000))
```

>[!IMPORTANT]
> It is necessary to  specify the column name present in the initial dataset. The implemented function relies on recognozinf these names, so it is recommended to avoid using names that may conflict with Python or Numpay reserved functions or keywords (e.g. log, sum, mean). Additionally, columns names containing special characters (such as +, -, *) should be avoided.
>To ensure proper functionlaity, it is recommended to use only letters, numbers and underscores, avoiding spaces in column names.

The name of the metric displayed in the graphs can be defined using the metric_name variable. For example:

```markdown
metric_name = 'Ranking Scale'
```

### Quick local use for non-experts

Perhaps MAYA is to run in Google Colaboratory, but if it is necessary you can download MAYA.ipynb and run in Jupythe with local resources. Also is possible clone the repository and follow these simple steps:

First we creat a new conda environment:
```markdown
# Here we name this conda environment 'MAYA' but you can call it anything you like
conda create --name MAYA --file environment.yml
```
We activate new conda environtment:
```markdown
conda activate MAYA
```
Now that we are in the environment we execute this comand
```markdown
# This is an example, with this intruction we are definig generate a chemical multiverse wit PCA ans t-SNE and just ECFP
./MAYA.py dataset='/content/example.csv', smiles_column_name='SMILES', target_activities=['Target_1', 'Target_2', 'Target_3'], MACCS=False, ECFP=True, MD=False, vPCA=True, t-SNE=True
``` 

[^1]:	Bertoni, M.; Duran-Frigola, M.; Badia-I-Mompel, P.; Pauls, E.; Orozco-Ruiz, M.; Guitart-Pla, O.; Alcalde, V.; Diaz, V. M.; Berenguer-Llergo, A.; Brun-Heath, I.; Villegas, N.; de Herreros, A. G.; Aloy, P. Bioactivity Descriptors for Uncharacterized Chemical Compounds. *Nat. Commun.* **2021**, *12* (1), 3932\.

