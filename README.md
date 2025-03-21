# Robust VEM for stress-assisted diffusion problems

This repository encapsulates all the experiments developed in this [pre-print](https://arxiv.org/abs/2401.09714).

The published version is available [here](https://epubs.siam.org/doi/10.1137/24M163640X)

## Description
In this work we present new results on a mixed nonlinear formulation for stress-assisted diffusion of a solute that interacts with an elastic material. Our study provides insights for a robust continuous analysis using parameter-dependent norms and the extended Babuˇska–Brezzi–Braess theory for perturbed saddle-point problems. The model describes the two-way coupling mechanism between the Herrmann formulation for linear elasticity and a mixed form of the reaction-diffusion equation consisting of diffusion-induced active stress and stress-dependent diffusion, respectively. The Banach fixed-point strategy with a small data assumption allows us to prove the well-posedness of the nonlinear coupled formulation. In addition, we propose stable VEM discretisations for the aforementioned non-linear problem. The a priori error analysis is studied and the theoretical robustness with respect to the physical parameters is confirmed numerically as well. Moreover, the applicability of our formulation is explored with a 2D model for the Lithiation of an anode, where our results coincide qualitatively with the literature.

## Getting Started

### Dependencies

MATLAB 9.14.0.2286388 (R2023a) Update 3

### Installing

Clone the repository with the following command using git:

$ git clone https://github.com/aerubianoma/vem_stress_assisted_diffusion

### Executing program

First of all, run the code init.m this allows you to set up MATLAB in the required folder. Then, run the desired experiment.

To change the manufacture solution or parameters go to [data](data) folder.

## Help

Contact the authors for bugs or recommendations.

## Authors

Rekha Khot

GitHub: [RekhaKhot](https://github.com/RekhaKhot).
Contact: rekha.khot@inria.fr

Andres E. Rubiano

GitHub: [aerubianoma](https://github.com/aerubianoma).
Contact: andres.rubianomartinez@monash.edu

Ricardo Ruiz-Baier

GitHub: [ruizbaier](https://github.com/ruizbaier).
Contact: ricardo.ruizbaier@monash.edu

## How to cite
Further details on the model and discretisation can be found in the following reference [[1]](https://arxiv.org/abs/2401.09714)
```
﻿@Article{Khot2025,
author={Khot, Rekha
and Rubiano, Andr{\'e}s E.
and Ruiz-Baier, Ricardo},
title={Robust Virtual Element Methods for Coupled Stress-Assisted Diffusion Problems},
journal={SIAM Journal on Scientific Computing},
year={2025},
month={Feb},
day={28},
publisher={Society for Industrial and Applied Mathematics},
volume={47},
number={1},
pages={A497-A526},
abstract={Abstract. This paper aims first to perform robust continuous analysis of a mixed nonlinear formulation for stress-assisted diffusion of a solute that interacts with an elastic material, and second to propose and analyze a virtual element formulation of the model problem. The two-way coupling mechanisms between the Herrmann formulation for linear elasticity and the reaction-diffusion equation (written in mixed form) consist of diffusion-induced active stress and stress-dependent diffusion. The two subproblems are analyzed using the extended Babu?ka?Brezzi?Braess theory for perturbed saddle-point problems. The well-posedness of the nonlinearly coupled system is established using a Banach fixed-point strategy under the smallness assumption on data. The virtual element formulations for the uncoupled subproblems are proven uniquely solvable by a fixed-point argument in conjunction with appropriate projection operators. We derive the a priori error estimates, and test the accuracy and performance of the proposed method through computational simulations. Reproducibility of computational results. This paper has been awarded the ?SIAM Reproducibility Badge: Code and data available? as a recognition that the authors have followed reproducibility principles valued by SISC and the scientific computing community. Code and data that allow readers to reproduce the results in this paper are available at https://github.com/aerubianoma/vem{\_}stress{\_}assisted{\_}diffusion and in the supplementary materials (vem{\_}stress{\_}assisted{\_}diffusion-main.zip [52.7MB]).},
issn={1064-8275},
doi={10.1137/24M163640X},
url={https://doi.org/10.1137/24M163640X}
}

```


