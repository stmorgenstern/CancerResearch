[![N|Solid|512x397](https://psor.uconn.edu/wp-content/uploads/sites/1972/2018/07/LabLogo_Graph_Full-768x238.png )](https://psor.uconn.edu/)
# Cancer Model Research
The cancer model research is the work of the cancer group within the PSOR laboratory at UCONN. This repository includes the most current version of the model implemented in Julia, the development model being used for active research in Spring 2021, and the older version of the model implemented in MATLAB.

| **Documentation**                                                |
|:-----------------------------------------------------------------:|
|[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/stmorgenstern/CancerResearch/blob/gh-pages/docslatest.md) | 


## Features and Capabilities
The current version of this model can be used to numerically model drug particles within the microvasculature and interstitium of a solid tumor based on the physics of the underlying transport phenomena. Specifically, our work provides functions for modeling velocity, pressure, and concentration profiles over time within a solid tumor. Furthermore, we have used these numerical models to formulate an optimization problem for effective permeability of drug particles.

## Work In Progress
- Formulate optimization problem using EAGO
- Create heat map of effective permeability in Julia
- Implement SciML surrogate model for concentration profile

# Contact
Samuel Degnan-Morgenstern, Undergraduate Chemical Engineering Researcher, samuel.morgenstern@uconn.edu

Dr. Matthew Stuber, Primary Investigator, matthew.stuber@uconn.edu
## References
1. Martin, J. D., Panagi, M., Wang, C., Khan, T. T., Martin, M. R., Voutouri, C., Toh, K., Papageorgis, P., Mpekris, F., Polydorou, C., Ishii, G., Takahashi, S., Gotohda, N., Suzuki, T., Wilhelm, M. E., Melo, V. A., Quader, S., Norimatsu, J., Lanning, R. M., Kojima, M., Stuber, M. D., Stylianopoulos, T., Kataoka, K., and Cabral, H. **Dexamethasone Increases Cisplatin-Loaded Nanocarrier Delivery and Efficacy in Metastatic Breast Cancer by Normalizing the Tumor Microenvironment.** *ACS Nano.* 13(6), 6396-6408 (2019).
