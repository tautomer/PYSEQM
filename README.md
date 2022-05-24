# [PYNEXMD: PYthon-based non-adiabatic molecular dynamics package](https://github.com/lanl/PYNEXMD)

[PYNEXMD](https://github.com/lanl/PYNEXMD) is a PYthon-based Nonadiabatic EXcited-State Molecular Dynamics package that leverages the [PyTorch](http://pytorch.org) and BML library for off-loading caculations on GPUs and/or heterogeneous exascale computing facilities. It provides build-in Semi-Empirical Quantum Mechanics method for electronic structures and interface to other first-principles packages (DFTB, PYSCF). It provides built-in interfaces for machine learning and efficient molecular dynamic engines with GPU supported. Several molecular Adiabatic and Non-adiabatic dynamics algorithms are implemented for facilitating dynamic simulations, inlcuding orginal and Extended Lagrangian Born-Oppenheimer Molecular Dynamics, geometric optimization, trajectory surface hopping, ab-initio multiple cloning, Ehrenfest dynamics.  

<hr/>

## Features:

* Interface with machine learning (ML) framework like [HIPNN](https://aip.scitation.org/doi/abs/10.1063/1.5011181) for ML applications and development.
* GPU-supported Molecular Dynamics Engine
* Off-load on heterogeneous exascale computing architecture via BML lib
* Efficient expansion algorithm [SP2](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.66.155115) for generating density matrix
* Stable and Efficient Extended Lagrangian Born Oppenheimer Molecular Dynamics ([XL-BOMD](https://aip.scitation.org/doi/full/10.1063/1.3148075))
* Trajectory surface hopping 
* Ab-initio multiple cloning
* Light-matter interaction
* Time-dependent density functional tight-binding 

## Installation:

```bash
git clone https://github.com/lanl/PYNEXMD.git (TO be relaesed there)
cd PYNEXMD
python setup.py install
```
or
```bash
pip install git+https://github.com/lanl/PYNEXMD.git
```

To enable GPU with CUDA, please refer to the Installation Guide on [PyTorch website](https://pytorch.org/)

## Prerequisites:
* PyTorch>=1.2

## optional package
* BML ()

## Usage:
see [```./doc/documentation.md```](./doc/documentation.md)

## Semi-Empirical Methods Implemented:
1. MNDO
2. AM1
3. PM3
4. DFTB

<hr/>

## Authors:

[Yu Zhang](mailto:zhy@lanl.gov), [Xinyang Li](mailto:lix@lanl.gov], [Guoqing Zhou](mailto:guoqingz@usc.edu), [Benjamin Nebgen](mailto:bnebgen@lanl.gov), Nicholas Lubbers, Anders M. N. Niklasson and Sergei Tretiak

## Citation:
TBA

## Acknowledgments:
Los Alamos National Lab (LANL), Center for Nonlinear Studies (CNLS), T-1

## Copyright Notice:

© (or copyright) 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.

## License:

This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
