# softpen_supplementary

This repo contains all materials to the paper:  *Maximum softly-penalized likelihood for mixed effects logistic regression*.

- **Data**: Contains the Culcita data and the conditional inference data, that is analysed in Sections 3 and 8 in the paper. 
- **Results**: Holds all numerical results and figures from analyses and simulations in the paper and the supplementary material.
- **Scripts**: All scripts to recreate the numerical results in the paper. 
	- **./Functions**: Contains the scripts `MSPAL.R` and `mv_MSPAL.R` which contain all functions used for the estimations and simulations as well as for the generation of figures and tables. 
	- `culcita_sep.R`: Generates the estimates in Table 2 of the main text and Table S1 of the supplementary material. 
	- `culcita_sim.R`: Conducts the simulation study in Section 3 of the main text and generates Figure 1.
	- `cond_inf_example.R`: Conducts all estimates and the simulation of Section 8 in the main text. Generates Figure 2 of the paper and Figure S1 of the supplementary material, as well as Tables 3 and 4 of the main text and table S2 in the supplementary material. 
	- `Sim_1.R`: Conducts the simulation study in Section S4.1 in the supplementary material and generates Figure S2 and Table S3 therein. 
	- `Sim_2.R`: Conducts the simulation study in Section S4.2 in the supplementary material and generates Figure S3 and Table S4 therein. 
- **softpen_supplementary.pdf**: Document with further material to the examples and simulations of the paper and additional simulation studies.

