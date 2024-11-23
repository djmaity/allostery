This repository contains the pipeline to create dataset of proteins with known allosteric and active sites

The code is in the notebooks directry and contains 4 Jupyter Python Notebooks. The Python packages used in the notebooks must be installed to succesfully run the notebooks. It is hihgly recommended to create a virtual environment with Anaconda or Miniconda. If not already intalled, install Anaconda or Miniconda following the instructions at https://docs.anaconda.com/anaconda/install/ and https://docs.anaconda.com/miniconda/miniconda-install/, respectively. Create the virtual environment with:

conda create -n allostery
conda activate allostery

conda install --yes --file requirements.txt

conda install pip
pip install -r requirements.txt

The nootbooks require an internet connection to automatically obtain the required data in many of its cells. The output is written to the output directory of the repository. The input and output locations can be customized in the notebooks. Please refer to the text cells in each notebook for further instructions. 

## 01_Clean_ASD_add_Active_Site_Residues.ipynb
Cleans the data from AlloSteric Database (ASD) and fetches the active site information from UniProt and Mechanism and Catalytic Site Atlas (M-CSA).
The final output of this notebook can be split into training, validation and test set for developing AI methods for allosteric site prediction.

Download the ASD_Release_202306_XF.tar.gz archive from https://mdl.shsmu.edu.cn/ASD/module/download/download.jsp?tabIndex=1 and extract it in the data directory of this repository before running the notbook so that the path is data/ASD_Release_202306_XF contain XML files

## 02_Compare_Training_Sets_of_Exisiting_Allosteric_Site_Prediction_Tools.ipynb
Identifies the UniProt and PDB IDs of the proteins in the training sets of existing allosteric site prediction tools.

## 03_Create_Common_Test_Set_for_Allosteric_Site_Prediction_Tools.ipynb
This notebook removes the proteins idenfitied by the notebook with prefix 02 from the proteins in the dataset obtained using notebook with predic 01.

Running the notebooks with prefixes 01, 02, and 03 in succession produces the test set without data leaks to benckmark exisinting allsoteric site prediction tools

## ASD_orthosteric_site.ipynb
Alternatively, the ASD provides a dataset of proteins with both allosteric and active site reisdues. However, this set is significantly smaller than the one obtained above. This notebook fetches the allosteric and orthosteric residues from ASD.