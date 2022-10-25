# ElectroPredictor
ElectroPredictor is a noncommercial Python applicaction designed for predicting the Mayr's electrophilicity index for a dataset of molecules. The script was implemented for Linux based operating systems. 
## Requisites
### Hardware
It is recommended to run ElectroPredictor using a device with 12 or more usable cores for an improved performance.
### Software
#### Python Libraries
ElectroPredictor was tested within a Anaconda Environment and it is highly recommended to run it in that way for avoiding
conflicts when importing the libraries.
  - OpenBabel and Pybel
  - Numpy
  - Pandas
  - python-weka-wrapper3
#### Applications
Both applications must be aded to the shell script of the operating system to freely execute them as commands in the terminal.
  - OpenJDK
  - MOPAC 22.0.4 
## How to set up ElectroPredictor?
1. Clone the repository files or download them as a zip file.
2. Extract ToMoCoMD QuBiLs-MIDAS Command Line Interface files in ToMoCoMD folder, but DO NOT replace headings.txt or the tomocomd_qubils.in. They are already configured to run ElectroPredictor in the most efficient way. DO NOT CHANGE THE FILES OR THE APPLICATION WILL NOT WORK. You can download the files for ToMoCoMD QuBiLs-MIDAS Command Line Interface (CLI) at: http://tomocomd.com/software/qubils-midas
3. Now it is ready to use!
## How to use ElectroPredictor?
1. Place an .sdf file in the main.py directory of ElectroPredictor. Make sure that the .sdf file is not empty.
2. Make sure all the libraries are updated and installed, and that you are able to run "mopac" and "java -jar" as commands in your Linux terminal.
3. Execute **main.py** 
4. The .csv with the results will be created at the **Results** folder.
## Additional information
The test sets used for the model validation are also supplied in the **TestSets** folder.
## References
- García-Jacas, C. R., Marrero-Ponce, Y., Vivas-Reyes, R., Suárez-Lezcano, J., Martinez-Rios, F., Terán, J. E., & Aguilera-Mendoza, L. (2020). Distributed and multicore QuBiLS-MIDAS software v2.0: Computing chiral, fuzzy, weighted and truncated geometrical molecular descriptors based on tensor algebra. Journal of Computational Chemistry, 41(12), 1209–1227. https://doi.org/10.1002/jcc.26167
- Nina Jeliazkova. (2017). ideaconsult/appdomain: ambit-appdomain v2.0.0 (ambit_appdomain-2.0.0). Zenodo. https://doi.org/10.5281/zenodo.265119
- Stewart, J. J. P., Klamt, A., Thiel, W., Danovich, D., Rocha, G. B., Gieseking, R. L., Moussa, J. E., Kurtz, H. A., Korambath, P., Merz, K. M., Jr., & Wang, B. (2022). MOPAC (Version 22.0.4) [Computer software]. https://doi.org/10.5281/zenodo.6511958


