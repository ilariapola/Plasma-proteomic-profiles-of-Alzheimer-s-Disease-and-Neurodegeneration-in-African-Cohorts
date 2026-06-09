# Sample R code from Nature Communications "Plasma proteomic profiles of Alzheimer's Disease and Neurodegeneration in African Cohorts"
This repository contains the R scripts used to reproduce the main analyses and figures presented 

Structure: 
### Demographics
Generation of demographic tables and descriptive statistics.

### Volcano_plots
Differential protein expression analyses:
- Amyloid-positive versus amyloid-negative participants.
- Cognitively impaired versus cognitively unimpaired participants within amyloid strata.

### Stacked histograms
Stacked bar plots showing the relative composition of functional protein categories among significantly altered proteins.

### Comorbidities analysis
Linear regression analyses assessing associations between grouped comorbidities and plasma protein levels

### LOESS
Protein trajectory visualizations across disease stages using z-scored protein concentrations relative to the CU- reference group.

## Software
R statistical software version 4.4.0.

Main packages:
- dplyr
- tidyr
- purrr
- broom
- ggplot2
- ggrepel
- gprofiler2
- table1

## License

MIT License.
