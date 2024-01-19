# Robust VEM for stress-assisted diffusion problems

This repository encapsulates all the experiments developed in this [publication](https://arxiv.org/abs/2401.09714).

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
@misc{khot2024robust,
      title={Robust virtual element methods for coupled stress-assisted diffusion problems}, 
      author={Rekha Khot and Andres E. Rubiano and Ricardo Ruiz-Baier},
      year={2024},
      eprint={2401.09714},
      archivePrefix={arXiv},
      primaryClass={math.NA}
}
```


