# Green's Matching

This README accompanies the paper "Green's Matching: An Efficient Approach to Parameter Estimation in Complex Dynamic Systems" by Jianbin Tan, Guoyu Zhang, Xueqin Wang, Hui Huang, and Fang Yao. The paper is accessible at [Journal of the Royal Statistical Society Series B: Statistical Methodology](https://academic.oup.com/jrsssb/advance-article-abstract/doi/10.1093/jrsssb/qkae031/7644665?redirectedFrom=fulltext&login=false).

## Citation
If you find our code or methodology useful, you may cite us as:

    @article{tan2024green,
      title={Green’s matching: an efficient approach to parameter estimation in complex dynamic systems},
      author={Tan, Jianbin and Zhang, Guoyu and Wang, Xueqin and Huang, Hui and Yao, Fang},
      journal={Journal of the Royal Statistical Society Series B: Statistical Methodology},
      pages={qkae031},
      year={2024},
      publisher={Oxford University Press UK}
    }

---
## Abstract
This repository contains the proofs, code, and a real dataset utilized in our study.

## 1. Data
### Abstract
The dataset employed in this study is sourced from the public domain and can be accessed at [FDA Software](https://www.psych.mcgill.ca/misc/fda/software.html). It comprises trajectories of Chinese script production, as detailed in sections 1.1.6 and 9.9 of "Dynamic Data Analysis" by Ramsay, James, and Giles Hooker (Springer, 2017).

### Availability
All necessary data to reproduce the results of this study are available.

### Data Dictionary
The real dataset, contained in the "dat.rda" file within the "Data" directory, includes two dynamic curves that represent the movement involved in producing Chinese script in both horizontal and vertical directions on a writing surface.

## 2. Code
### Abstract
Our study employs Green's matching technique for parameter estimation in dynamic systems, compared to other methods such as:
- [Gradient matching and Integral Matching](https://link.springer.com/content/pdf/10.1007/978-1-4939-7190-9.pdf)
- [Generalized Smoothing Approaches](https://academic.oup.com/jrsssb/article/69/5/741/7109525)
- [Manifold-Constrained Gaussian Processes](https://www.pnas.org/doi/abs/10.1073/pnas.2020397118)

### Reproducibility
- "Comp_sim.R" was used to generate the results presented in Table 1 of the main text.
- "Plot_sim.R" was executed to illustrate estimation biases and variances of the model parameters shown in Figure 2.
- Curve reconstructions and confidence interval calculations for Figure 3 were conducted using "Curve_inf.R" and "Conf_inf.R". Adjustments to the "rat" parameter were made to produce results shown in Figures 2-5 and Table 1 of the Supplementary Materials.
- The equation discovery process described in Section 5 was implemented in "Running.R" within the "Equation_discovery" directory.
- Additional simulation studies mentioned in Section 4.2 and data illustrations in the Supplementary Materials were carried out using "Comp_sim_hand.R" and "Data_illustration.R", respectively.

## 3. Supplementary Materials
The proofs supporting the methodologies used in this study are detailed in "SM.pdf".

