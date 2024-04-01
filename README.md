# Green's matching: an efficient approach to parameter estimation in complex dynamic systems

This file contains the data and code used in the paper "Green's matching: an efficient approach to parameter estimation in complex dynamic systems" by Jianbin Tan, Guoyu Zhang, Xueqin Wang, Hui Huang, and Fang Yao.

# Abstract
This file contains the proof, codes and a real dataset used in this article. 

## 1. Data
### 1) Abstract

The data underlying this article are derived from a source in the public domain, available at https://www.psych.mcgill.ca/misc/fda/software.html. The dataset contains trajectories in producing a Chinese script. Additional information on the Chinese handwriting data can be found in Sections 1.1.6 and 9.9 of "Ramsay, James, and Giles Hooker. Dynamic data analysis. Springer, 2017."

### 2) Availability
The data to reproduce our results are available.

### 3) Data dictionary
The real dataset is included in the "dat.rda" file within the "Data" directory. It consists of two dynamic curves that capture the movement of producing a Chinese script in both horizontal and vertical directions on a writing surface.

----
## 2. Code
### 1) Abstract
We employed Green's matching for parameter estimations in dynamic systems compared with other methods, including the gradient matching of different orders, the generalized smoothing approach, and the manifold-constrained Gaussian processes.

### 2) Reproducibility
- The results in Table 1 (in the main text) were produced by running "Comp_sim.R".
- The illustration of estimation biases and variances of the model parameters in Figure 2 was conducted by running "Plot_sim.R".
- The curves' reconstructions and confidence intervals calculations in Figure 3 were conducted by running "Curve_inf.R" and "Conf_inf.R", respectively, in which we changed the parameter "rat" into other values to obtain the results in Figures 2-5 and Table 1 (in Supplementary Materials).
- The equation discovery in Section 5 was implemented by "Running.R" in the "Equation_discovery" directory.
- The additional simulation studies in Section 4.2 and the data illustration in Supplementary Materials were obtained by running "Comp_sim_hand.R" and "Data_illustration.R", respectively.

----
## 3. Supplementary materials
The proof of this article is included in the "SM.pdf".
