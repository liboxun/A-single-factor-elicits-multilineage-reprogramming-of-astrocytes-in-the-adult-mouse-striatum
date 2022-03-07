# A-single-factor-elicits-multilineage-reprogramming-of-astrocytes-in-the-adult-mouse-striatum

This contains all the code needed to reproduce the scRNA-seq analyses and figures in the publication, 'A single factor elicits multilineage reprogramming of astrocytes in the adult mouse striatum', PNAS, 2022, Vol. 119  No. 11 e2107339119 (https://doi.org/10.1073/pnas.2107339119). 

`0.Run_Scrublet.ipynb` is used to process the gene-cell expression matrices as a result of the cell-ranger alignment pipeline. 

`1.MainAnalysis_AllCells.ipynb` is used to make parts of Fig. 6 and S18. 

`2.MainAnalysis_ReprogrammingAlone.ipynb` is used to make parts of Fig. 5 and S13. 

Code for all Gene Ontology-related analysis is in `/GOEA`.

Code for running pySCENIC (Fig. 5G) is in `/pySCENIC`.

Code for comparing DLX2-induced reprogramming and adult neurogenesis (Fig. S19) is in `/pySCENIC`. 

Everything else is in `/Misc`. 

For any questions or problems, open an issue or contact the author at boxun.li@duke.edu
