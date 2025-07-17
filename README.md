## SLIDE

![Graphical Abstract](https://github.com/user-attachments/assets/c41a1f86-8e47-45f0-b4d6-71cf0cea32a3)

# Introduction to SLIDE

Sliding identification by Local Integration of Displacements across Edges (**SLIDE**) is a methodology to locally identify and quantify Grain Boundary Sliding (GBS) and Grain Boundary Opening (GBO). Calculation of the GBS/GBO magnitude is performed through the matching of Digital Image Correlation (DIC) displacement gradient fields to the local kinematics of the grain boundary, obtained from EBSD. In addition, GBO is identified using In-Beam SE images, and the 3D sliding magnitude is obtained through Optical Profilometry fields. More details can be found in [**this preprint**](https://arxiv.org/abs/2504.04898). The SLIDE framework works in similar fashion to the Slip Systems based Identification of Local Plasticity (SSLIP) method ([**see this paper**](https://www.sciencedirect.com/science/article/pii/S1359645422008795?via%3Dihub)); the two methods can be used in conjunction for combined identification of intragranular crystallographic slip and grain boundary deformation.

The **SLIDE** function library is written in [**MATLAB**](https://mathworks.com/products/matlab.html) and uses several functionalities of the MATLAB-based crystallographic toolbox [**MTEX**](https://mtex-toolbox.github.io). 

**WARNING: It is recommended to use Mtex 5.11.2.**

It is important to use aligned EBSD/DIC data. See the following repository for an alignment framework: [**NanoMech_Alignment_Matlab**](https://github.com/Tijmenvermeij/NanoMech_Alignment_Matlab).

Please report any bugs you encounter.

# Authors
**SLIDE** has been conceptualized and created by [**Casper Mornout**](https://research.tue.nl/en/persons/casper-mornout), **Gert-Jan Slokker**, **Dennis König**, **Tijmen Vermeij** and **Johan Hoefnagels**.

# How to cite SLIDE
If you have applied the SLIDE code to your research, please cite this open-access paper as your reference:
C.J.A. Mornout, G. Slokker, T. Vermeij et al., SLIDE: Automated identification and quantification of grain boundary sliding and opening in 3D, Scripta Materialia, https://doi.org/10.1016/j.scriptamat.2025.116861
