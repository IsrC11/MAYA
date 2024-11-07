# MAYA
<div align='center'>
  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/17SSd2BuBfMffRJKJwfYrDPwK3khvvAIj?usp=sharing)


<p align='justify'>
Multiple ActivitY Analyzer (MAYA) is designed to automatically construct a chemical multiverse, generating multiple visualizations of chemical spaces described by structural descriptors such as MACCS keys (166 bits), ECFP 4 and 6, molecular descriptors with pharmaceutical relevance as well as implementing biological descriptors. These representations are integrated with various visualization techniques for automated analysis, focusing on the analysis of structure - multiple activiy/property relationships (SMARTs).
  
<p align='center'>
<img src="https://github.com/ApSirius/Autimated-Analysis-of-Structure-Multiple-Property-Relationships/blob/d35816f1e1c1e1da9aa98790b7fc91065f5f6162/Chemical%20multiverse.png" alt="Process" width="600">
<p align='justify'>
</div>
<p align='justify'>
MAYA has been developing as a user-friendly, open-source tool that automates the construction of chemical spaces by integrating diverse molecular representations to provide a more comprehensive description of the structural, chemical, and functional characteristics of a given set of molecules described by their SMILES notation and an associated activity/property, and the tool supports various file formats (CSV, TSV, XLSX, JSON and XML), requiring only the specification of a few parameters related to the database in use and the desired representations. Additionally, MAYA includes options for customizing the visualizations.
<p align='justify'>
The generated visualizations are interactive, enhancing the understanding of the displayed data. They offer a 2D view of the molecular structure, along with the obtained variability values from PCA and their SMILES notation. Customization features allow users to adjust the size, shape, and transparency of data points, as well as the ability to modify the color palette.

<p align='justify'>
The script consist in a funtion that automatically implement:
  
1. Data curation
2. Descriptors calculation
3. Tanimoto simmilarity calculation
4. Dimensionality reduction
5. 2D interactive visualization

### How use MAYA?

>[!IMPORTANT]
>It is important to be sure our dataset contains the follow information:
>1. Smiles notation
>2. Identificator
>3. Value of activity or properties
>Depending of the interest of the user it is possible select the descriptors and dimensionality reduction thecniques to use. Defining the variables as True or False is possible disable their calculation.

### Example of usage
```markdown
# This is an example
chemical_multiverse(dataset='/content/example.csv', smiles_column_name='SMILES', target_activities=['Target_1', 'Target_2', 'Target_3'], MACCS=Falce, ECFP=True, MD=Falce, vPCA=True, t-SNE=True )
```
See this [notebook](https://github.com/IsrC11/MAYA/blob/d2ca032691cf98c1ae805c567a6b4508bf5dc168/Local_usage.ipynb) for more detailed usage

### Why use MAYA?
<p align='justify'>
To perform an automated analysis of your database annotated with any activity, property, or score by constructing a chemical multiverse focused on a deeper understanding of multiple structure-activity relationships. 
<p align='justify'>
You can customize the descriptors and techniques used depending on the required focus. You can select which descriptors you want to use, and you can also input a similarity matrix of any desired descriptor, allowing its integration into the generated visualizations.
<p align='justify'>
Access to well-documented code is provided, covering database curation processes, similarity calculations, and dimensionality reduction techniques.

### Usage
1. Google Colaboratory <br> The easiest way to use the script is ti open it in [Google Colaboratory](https://github.com/IsrC11/MAYA/blob/d2ca032691cf98c1ae805c567a6b4508bf5dc168/Local_usage.ipynb). The only thing needed is a Google account.
2. Local installation <br> You can also setup your own local environment if you do not want to run the script through a Google service.

### Additional Information
<p align='justify'>
MAYA current supports Pythob 3.10

rdkit (2022.09.05)

matplotlib (3.7.1)

pandas (2.1.4)

seaborn (0.13.1)

sklearn (1.3.2)

### Funding
Research contained in this package was supported by the Consejo Nacional de Humanidades, Ciencia y Tecnología (CONAHCYT) for the scholarship No. CVU 1340927
