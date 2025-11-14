# MAYA
<div align='center'>
  
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/drive/1fsHYMMJ5CgaC9NuVRN4E0zVC9qz28XQ0?usp=sharing)


<p align='justify'>
Multiple ActivitY Analyzer (MAYA) is an open-source Python tool designed to automate the generation of chemical multiverses for comprehensive analysis of structure-activity/property relationships/associations. MAYA integrates multiple molecular representations, including structural descriptors (e.g., MACCS keys, ECFP4, ECFP6, and MAP4), physicochemical molecualr descriptrs, and biological descriptors, to construct diverse chemical spaces. These spaces are visualized through interactive 2D plots, enabling deeper insights into molecular characteristics and their relationships with activities or properties

<p align='center'>
<img src="https://github.com/IsrC11/MAYA/blob/324aefd424da8d25b9dcc5ef2caba67743437a09/Wflow.png" alt="Process" width="1000">
<p align='justify'>
</div>
<p align='justify'>
MAYA is user-friendly and supports various input file formats (CSV, TSV, XLSX, JASON, and XML), requering only a dataset with SMILES notation, molecualr identifiers, and associated activity/property values. Users can customize the analysis by selecting specific descriptors and demensionality reduction techniques (e.g., PCA) through simple parameters setting. The tool automates the following key processes: 
<p align='justify'>
The script consist in a funtion that automatically implement:
  
1. Data curation: Ensures high-quality data by validating and preprocessing molecualr datasets.
2. Descriptors calculation: Computes structural, physicochemical, and biological descriptors for molecular characterization.
3. Tanimoto simmilarity calculation: Quantifies molecular similarity to support chemical space analysis
4. Dimensionality reduction: Applies techniques like PCA and t-SNE to reduce complexity while preserving meaningful patters.
5. 2D interactive visualization: Generates customizable, interactive plots displaying molecualr structures, PCA-derived variability, SMILES notaton, and user-defined visual attributes (e.g., point size, shape, color palette, and transparency)

MAYA's visualizations are designed to enhance intepretability, offering researchers a clear and interactive way to explore chemical spaces. The tools is particularly suited for chemoinformatics applications, such as drug discovery, where understanding complex structured-activity relationships is critical.

[Here](https://github.com/IsrC11/MAYA/blob/cc7d8d7f59947269491ba37b3f772eeca4d81741/User_Guide.md) you can find more detailed information about how MAYA works 

### How use MAYA?

>[!IMPORTANT]
>It is essential to ensure our dataset contains the following information:
>1. Smiles notation: Molecular representation
>2. Identifier: Unique molecule identifiers
>3. Activity or property values: Quantitative or categorical data for analysis
>
>Users can customize the analysis by enabling or disabling specific descriptors and dimensionality reduction techniquese via Boolean parameters (e.g., True or False) Detailed instructions and examples are available in the documentation.

#### Example of usage
```markdown
# This is an example
from maya_chem import MayaAnalyzer, MayaConfig
import numpy as np

# Crear configuración
config = MayaConfig(data_path="TEST.xlsx")

# Actualizar claves relevantes
config.data.update({
    "id_col": "molregno",
    "smiles_col": "canonical_smiles",
    "activities": ["standard_value"]
})
config.analysis["fingerprint"] = ["morgan", 'maccs']
config.analysis["reduction_method"] = ['pca']
config.viz["output_dir"] = "/content/MAYA/colab_results"

# Correr pipeline
analyzer = MayaAnalyzer(config)
results = analyzer.run(color_by='standard_value')
figs = analyzer.visualize(interactive_mode=True, color_by='standard_value')

print("✅ Pipeline completed. Check '/content/MAYA/colab_results' for outputs.")

```

<p align='center'>
<img src="https://github.com/IsrC11/MAYA/blob/44661380db69df87131c8dd5315a24debe5e6020/Example.jpg" alt="Process" width="1000">
<p align='justify'>
</div>


See this [notebook](https://github.com/IsrC11/MAYA/blob/cc7d8d7f59947269491ba37b3f772eeca4d81741/User_Guide.md) for more detailed usage

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
