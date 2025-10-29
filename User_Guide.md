# **User Guide**

<p align='justify'>
MAYA is a tool designed to automate the construction and visualization of a chemical multiverse from a set of compounds. It enables multiple similarity calculations, dimensionality reduction, quality metrics, and molecule-on-hover visualizations based on chemical descriptors and biological activities.

*Key Features*

Feature     |    Descriptions
 :---:      |     :---:
SMILES Preprocessing  |  Normalization, validation, and canonicalization
Fingerprints          |  MACCS Keys (166 bits), Morgan/ECFP (1024, 2048), MAP4
Similarity            |  Parallel Tanimoto matrix (via Joblib)
Dimensionality Reduction  |  PCA, t-SNE, UMAP
Physicochemical Properties | MolWt, LogP, HBA, HBD, TPSA
Quality Metrics       |  Trustworthiness, Distancce Correlation
Vsualization          |  Static (PNG) + Interactive (Dash + Molplotly)
CLI & API             |  Run from terminal or Python

---

### Dataset Format  
<p align='justify'>
 
>[!NOTE]
>The dataset must be annotated with SMILES notation (required), a unique idenntifier ID (required), and at least one activity or property (optional). Supported formats include: CSV, XLSX, TSV.

Comp_ID   |    Smiles    |    Activity_1
 :---:      |     :---:      |      :---:
Comp_137  |  NC(=O)c1cccc(c1)c2cnc3[nH]cc(c4ccccc4)c3c2  |  0.001


### Analysis Configuration

- **data_path** (str): Input file path

- **output_dir** (str): Output folder. Default: ./results

- **fingerprint** (str/list): List of fingerprints for compute. Default: ['maccs']

- **reduction_method** (str/list): Visualization technique to use. Default: ['pca']

- **color_by** (str): Column to color points. Default: 'LogP'

- **palette** (str): Plotly/Matplotlib palette. Default: 'RdBu_r'

  #### Customize columns
```markdown
config.data['id_col'] = 'Molecule_ID'
config.data['smiles_col'] = 'Structure'
config.data['activities'] = ['pIC50_A', 'pIC50_B']
```
 
 ## Key parameters


 #### Fingerprints

- **Morgan (ECFP)**: 1024/2048 bits (Recommended)

- **MACCS** (float): 166 bits. Use for fast screening

- **MAP4** (str)=Variable to define the property to be vis

```markdown
# In descriptors.py compute_morgan_fingerprint(radius=3)
```

#### Dimensionality Reduction

Method   |    Key Parameters    |    Recommended Range   |    Effect
 :---:      |     :---:      |      :---:    |      :---:
PCA         |  n_compontes  |  2-3  |   Linear, Interpretable
t-SNE       |  perplexity   |  5-50  |  Local structure
UMAP        |  n_neighbors  |  5-200 |  Local + Global balance

How many neighbors does each point consider?

Value   |    Effect
 :---:      |     :---:
5-10         |  Small, tight clusters
30       |  Balacned (default)
50-100        |  Global trends

##### Customize
```markdown
from sklearn.mainfold import TSNE
tsne = TSNE(perplexity=50, n_iter=5000)
```

### Physicochemcial Properties

Value   |    Meaning
 :---:     |     :---:
MolWt      |  Molecular Weight
LogP       |  lipophilicity
HBA        |  Hydrogen bond acceptors
HBD        |  Hydrogen bond donors
TPSA       |  Polar Surface Area


---

### Best Practices

Topic   |    Reocommendation
 :---:     |     :---:
SMILES      |  Always canonicalize
Size       |  >50k compounds
Reproducibility        |  random_state=42
Interpretation        |  PCA (global) + t-SNE/UMAP (local)
TPSA       |  Polar Surface Area

---
>[!IMPORTANT]
>1. Data visualization: Ensure the dataset does not contain null values in key columns (ID and smiles_column_name).
>2. Parallelization: Use n_jobs=/1 to utilize all available cores and speed calculations.
>3. Metric Selection: Carefully defined the evaluation metric to obtain meaninful results.
>4. Visualization: Adjust ˋsize_point_paletteˋ, and ˋpoint_shapeˋ to costumize plots as needed.

---

## Relevance of some variables
<p align='justify'>
   
##### **radius**  
<p align='justify'>

Defines the radius of the circle centered on each atom within a molecule. Typically this radius corresponds to 2 or 3 bonds, allowing the mapping of the immediate chemical environment of each atom. This variable is essential for capturing local chemical interactions that contribute to molecular activity.

#### **signaturizer\_code**  
A list containing the names of the different bioactive descriptors[^1] used. These descriptors can be expanded or reduced depending on the analysis requirements. The included descriptors are:

Level     |   1  |   2  |  3  |   4  |   5
  :---:     |  :---:  |  :---:  |  :---: |  :---:  |  :---:
Chemistry (A) |  2D Fingerprints  |  3D Fingerprints  |  Scaffolds |  Structural keys  |  Physicochemistry
Targets (B)   |  Mechanisms of action  |  Metabolic genes  | Crystals |  Binding  |  HTS bioassays
Networks (C)  | Small molecule roles | Small molecule pathways | Signaling pathways | Biological proceses | Interactome
Cells (D)     |  Transcription  |  Cancer cell lines  | Chemical genetics | Morphology | Cell bioasssays
Clinics (E)   | Therapeutic areas |  Indications  | Side effects |  Diseases & Toxicology |  Drug-Drug Interactions

#### **palette**  
This variable is used to modify the color palette applied in the graphs. It supports two main configurations\>

1. Predefined Matplotlib palettes: Users can select one of the color palettes already available in the Matplotlib library, such as ‘viridis’, ‘plasma’ or ‘coolwarm’.  
2. Custom list: Alternatively, a custom list of colors can be defined using RGB or HEX codes (e.g. \#FF5733 for orange).

This flexibility allows users to tailor visualizations to specific needs or aesthetic requirements.

#### **perplexity**  
The perplexity value is closely related to the resulting visualizations. This parameter represents an approximation of the number of nearest neighbors for each point, balancing the weight between local and global features in the analysis. It is recommended to set a perplexity value between 5 and 50 depending on the database size, ensuring that the value is always less than the total size of the database

## Color Scale
<p align='justify'>
The visual representation of a property or activity is implemented to facilitate the identification of relevant patterns, differences, or trends. Additionally, weighting and penalizations considerations can be implemented to enhance the analysis copabilities based on user requeriments. This allows the generation of a ranking or evaluation matric based on the activities or properties available in the database.
  

### Quick local use for non-experts

Perhaps MAYA is to run in Google Colaboratory, but if it is necessary you can download MAYA.ipynb and run in Jupythe with local resources. Also is possible clone the repository and follow these simple steps:

First we creat a new conda environment:
```markdown
# Here we name this conda environment 'MAYA' but you can call it anything you like
conda create --name MAYA
```
We activate new conda environtment:
```markdown
conda activate MAYA
pip install -r requirements.txt
```
Now that we are in the environment we execute this comand
```markdown
# This is an example, with this intruction we are definig generate a chemical multiverse wit PCA ans t-SNE and just ECFP
./MAYA.py dataset='/content/example.csv', smiles_column_name='SMILES', target_activities=['Target_1', 'Target_2', 'Target_3'], MACCS=False, ECFP=True, MD=False, vPCA=True, t-SNE=True
``` 

[^1]:	Bertoni, M.; Duran-Frigola, M.; Badia-I-Mompel, P.; Pauls, E.; Orozco-Ruiz, M.; Guitart-Pla, O.; Alcalde, V.; Diaz, V. M.; Berenguer-Llergo, A.; Brun-Heath, I.; Villegas, N.; de Herreros, A. G.; Aloy, P. Bioactivity Descriptors for Uncharacterized Chemical Compounds. *Nat. Commun.* **2021**, *12* (1), 3932\.

