# scrnaseq
Molofsky lab database of single-cell rnaseq workflow (https://www.annamolofskylab.org/)




Welcome! Check out the "Vignettes" folder for worked examples of each R Script (scripts can be found in the "Rmd" folder). Start with the A-B-C-D scripts, the remainder can be used as needed to generate different kinds of plots, visible in the associated vignette. 




You will need to edit the first 1-3 chunks of code for each script to be specific to your project and your directory/folder system. Feel free to email me (leah.dorman - at - ucsf.edu) if you have any problems loading your data or getting a script to run properly. 






Troubleshooting tips: 
You may have to install new packages the first time - usually BiocManager::install("package-name") will work. 
If a previously installed package command is not working, run library(package-name). 
If your data will not load, the path (directory) is likely incorrect. Try right-clicking on the file itself on your computer and hitting the "option" key to "copy as Pathname" and directly paste that into the R script for any file that won't load properly. 
You can always subset your data and re-run the same script with the new data - just change the  "name" in the top chunk so that any saved spreadsheets or plots don't overwrite the old ones. You can of course put the old ones in a different folder to be safe. 






Just a note - "Cellphonedb" and the "spliced gene/ratio" work with spreadsheets you will have to create from your data using Python code for cellphonedb and scvelo, respectively. These scripts are in the "Python" folder but may not be fully generalizable yet and take slightly more editing to work with your system. 

Thanks to: 
Satija Lab (Seurat package) https://satijalab.org/seurat/;
Matrisome Project http://matrisomeproject.mit.edu/;
Gabriel Hoffman (Variance Partition) https://www.bioconductor.org/packages/release/bioc/html/variancePartition.html;
Teich Lab (CellPhoneDB) https://www.cellphonedb.org/;
Kevin Blighe (Enhanced Volcano) http://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html;
UC Davis Bioinformatics Core;
Barres Lab/Ye Zhang (Bulk cortical RNA-Seq) http://brainrnaseq.org;
Tirosh/Izar/Regev/Garraway, Science 2016 (cell cycle assignment) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4944528/;

