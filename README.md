# MAYA
<div align='center'>
  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/17SSd2BuBfMffRJKJwfYrDPwK3khvvAIj?usp=sharing)
<p align='center'>

**MAYA** (Multiple Activity Analyzer) is designed to automatically construct a chemical multiverse, generating multiple visualizations of chemical spaces described by structural descriptors such as MACCS keys (166 bits), ECFP 4 and 6, and molecular descriptors with pharmaceutical relevance, as well as implementing biological descriptors. These representations are integrated with various visualization techniques for automated analysis, focusing on the analysis of structure - multiple activiy/property relationships.

![Process](https://github.com/ApSirius/Autimated-Analysis-of-Structure-Multiple-Property-Relationships/blob/d35816f1e1c1e1da9aa98790b7fc91065f5f6162/Chemical%20multiverse.png)
---
</div>
<p align='justify'>
MAYA is developed as a user-friendly, open-source tool that automates the construction of chemical spaces by integrating different representations to provide a more comprehensive description of the structural, chemical, and functional characteristics of a set of molecules described by their SMILES notation and an associated activity/property, supporting various file types (CSV, TSV, XLSX, JSON and XML), requiring only the specification of a few parameters related to the database in use and the desired representations. Additionally, it includes options for customizing the visualizations.
<p align='justify'>
The generated visualizations are interactive, allowing for a better understanding of the displayed data. They provide a 2D view of the structure, as well as the obtained variability values and their SMILES notation. Customization features are included, enabling the modification of the data's size, shape, and transparency, as well as the ability to change the color palette.

<p align='justify'>
The script consist in a funtion that automatically implement:
  
1. Data curation
2. Descriptors calculation
3. Tanimoto simmilarity calculation
4. Dimensionality reduction
5. 2D interactive visualization

The user only need to asign some variables related to the dataset:

1. Name of the ID column
2. Name of the SMILE column
3. Name of Activity/Property columns
4. Define the descriptors and kind of visualization are required

### Why use MAYA?
<p align='justify'>
To perform an automated analysis of your database annotated with any activity, property, or score by constructing a chemical multiverse focused on a deeper understanding of multiple structure-activity relationships. 
<p align='justify'>
You can customize the descriptors and techniques used depending on the required focus. You can select which descriptors you want to use, and you can also input a similarity matrix of any desired descriptor, allowing its integration into the generated visualizations.
<p align='justify'>
Access to well-documented code is provided, covering database curation processes, similarity calculations, and dimensionality reduction techniques.

### Usage
1. Google Colaboratory <br> The easiest way to use the script is ti open it in Google Colab. The only thing needed is a Google account.
2. Local installation <br> You can also setup your own loclal Python environment if you do not want to run the script through a Google service.

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
