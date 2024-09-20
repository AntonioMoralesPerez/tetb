# Transversality-Enforced Tight-Binding (TETB) model

**Implemented by Antonio Morales-Pérez, Chiara Devescovi, Yoonseok Hwang, Mikel García-Díez, 
    Barry Bradlyn, Juan L. Mañes, Maia G. Vergniory, and Aitzol García-Etxarri**

In this package we define the essential functions necessary to implement the *enumeration algorithm*, 
which returns all possible decompositions in terms of elementary bands representations (EBRs) 
of a "symmetry vector" of a bundle of isolated photonic bands, given a certain number of auxiliary bands.

The analytical framework underlying our implementation of the TETB model is detailed in Ref. [1]. 
We will use the notation for labelling the space groups (SGs), irreducible representations (irreps) 
and EBRs presented in Ref. [3].

## Installation

Clone the repository and install the package by navigating to the working directory and running

````
pip install .
````

## Requirements

It requieres three additional packages:

1. [NumPy](https://github.com/numpy/numpy.git)
2. [SymPy](https://github.com/sympy/sympy.git) (v1.12)
3. [OR-Tools - Google Optimization Tools](https://github.com/google/or-tools.git) (v9.10.4067)

## Examples

We provide a set of Jupyter notebooks in `./examples/` for the following cases:

1. The specific examples studied in Ref. [1] under space groups (SGs) #224 and #221
2. All first minimal solutions of SG #224, #221, #99 and #2 presented in Ref. [2]

As a first step, it is necessary to define by hand all the EBRs and irreps of the specific 
SG of the system under study. To accomplish this we forward the user to the tables shown in 
Ref. [3-5]. There, one can find all EBRs and all irreps tabulated for all magnetic and non-magnetic 
SGs. More details are given within the examples.

## References

Please cite the following paper when using this package:

[1] A. Morales-Pérez, C. Devescovi, Y. Hwang, et al., “Transversality-enforced tight-binding 
model for 3d photonic crystals aided by topological quantum chemistry,” 
arXiv preprint arXiv:2305.18257 (2023).

[2] T. Christensen, H. C. Po, J. D. Joannopoulos, and M. Soljačić, “Location and topology of 
the fundamental gap in photonic crystals,” Phys. Rev. X 12, 021066 (2022).

[3] B. Bradlyn, L. Elcoro, J. Cano, M. G. Vergniory, Z. Wang, C. Felser, M. I. Aroyo & B. A. 
Bernevig
"Topological quantum chemistry"
Nature 547, 298-305 (2017).

[4] [Bilbao Crystallographic Server -> BANDREP](https://www.cryst.ehu.es/cgi-bin/cryst/programs/bandrep.pl)

[5] [Bilbao Crystallographic Server -> Representations SG](https://www.cryst.ehu.es/cgi-bin/cryst/programs/representations.pl?tipogrupo=spg)

## Acknowledgments

A.G.E., A.M.P, M.G.D., M.G.V. and C.D. acknowledge support from the Spanish Ministerio de 
Ciencia e Innovación (PID2022-142008NB-I00).
A.G.E., A.M.P, and C.D. also acknowledge support from the Basque Government Elkartek program 
(KK- 533 2023/00016) and the Gipuzkoa Provincial Council within the QUAN-000021-01 project.
A.G.E. and M.G.V. acknowledge funding from the IKUR Strategy under the collaboration agreement 
between Ikerbasque Foundation and DIPC on behalf of the Department of Education of the Basque 
Government, Programa de Ayuda de Apoyo a los agentes de la Red Vasca de Ciencia, Tecnolog\'ia e 
Innovaci\'on acreditados en la categor\'ia de Centros de Investigaci\'on B\'asica y de Excelencia 
(Programa BERC) from the Departamento de Universidades e Investigaci\'on del Gobierno Vasco and 
Centros Severo Ochoa AEI/CEX2018-000867-S from the Spanish Ministerio de Ciencia e Innovaci\'on.
M.G.V. thanks support to the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) 
GA 3314/1-1 – FOR 5249 (QUAST) and partial support from European Research Council (ERC) grant 
agreement no. 101020833. The work of JLM has been partly supported by the Basque Government 
Grant No. IT1628-22 and the PID2021-123703NB-C21 grant funded by MCIN/AEI/10.13039/501100011033/ 
and ERDF; ``A way of making Europe”. The work of B.B. and Y.~H. is supported by the Air Force 
Office of Scientific Research under award number FA9550-21-1-0131. Y.~H. received additional 
support from the US Office of Naval Research (ONR) Multidisciplinary University Research
Initiative (MURI) grant N00014-20-1-2325 on Robust
Photonic Materials with High-Order Topological Protection. C.D. acknowledges financial support 
from the MICIU through the FPI Ph.D. Fellowship CEX2018-000867-S-19-1. M.G.D. acknowledges 
financial support from the Government of the Basque Country through the pre-doctoral fellowship 
PRE\_2022\_2\_0044.

## License

This project is released under the [GNU GENERAL PUBLIC LICENSE](./LICENSE)