# mNEPv - Monotone Nonlinear Eigenvector Problems

Dated 		04-20-2023

This folder contains the MATLAB implementation of the algorithms 
used in the research paper 

>*Variational Characterization of Monotone Nonlinear Eigenvector Problems and Geometry of Self-consistent Field Iteration*
by Zhaojun Bai and Ding Lu, 2022.
(Manuscript available at: https://doi.org/10.48550/arXiv.2211.05179 )


## Description

A Monotone Nonlinear Eigenvector problems (mNEPv) is defined as 
$$H(x) x = \lambda x,$$
where $H(x)$ is a Hermitian matrix-valued function of the form 
$$H(x):= \sum_{i=1}^m h_i(x^HA_ix) A_i,$$
and $A_1,\dots,A_m$ are $n$-by-$n$ Hermitian matrices, 
$h_1,\dots,h_m$ are differentiable and non-decreasing functions over $\mathbb R$. 
The goal is to find a unit-length vector $x\in\mathbb C^n$ and a scalar $\lambda\in\mathbb R$ 
satisfying $H(x)x=\lambda x$ and, furthermore, $\lambda (= x^H H(x) x)$ is the largest eigenvalue of $H(x)$.

The algorithm implemented is the *self-consistent-field iteration* (SCF) with
local acceleration, which ensures global convergence to a solution, as described in the paper.
This package contains all the example data and routines used in the study to facilitate its reproducibility.


## Contents

- mNEPv solver
	- scf.m:			SCF iteration with and without local acceleration

- Examples 
	- numrd2d.m:  		mNEPv from numerical radius computation  
	- dhdae1.m: 		mNEPv from distance problems of linear dHDAE systems
	- dhdae2.m:			mNEPv from distance problems of quadratic dHDAE systems
	- tensorappr.m:  	mNEPv from Tensor rank-1 approximation

- Other files: numrange.m (plot numerical range), goescf:.m (geometry of SCF), and data files.


## External

The MATLAB toolbox Manopt is used in dhdae2 for comparison. The
comparison is deactivated by default. To use it please first download
and install the [Manopt](https://www.manopt.org/).


## Contact 

For questions, please contact Ding.Lu@uky.edu  

