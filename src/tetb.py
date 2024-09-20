# necessary additional packages for the code to work
import numpy as np
import sympy as sp
from ortools.sat.python import cp_model

fixedGamma2 = (
    [1]
    + list(np.arange(3, 5 + 1, 1))
    + list(np.arange(6, 9 + 1, 1))
    + list(np.arange(25, 46 + 1, 1))
    + list(np.arange(75, 80 + 1, 1))
    + list(np.arange(143, 146 + 1, 1))
    + list(np.arange(156, 161 + 1, 1))
    + list(np.arange(168, 173 + 1, 1))
    + list(np.arange(183, 186 + 1, 1))
)  # SGs with fixed surrogate representation with \Gamma_2 as trivial irrep

fixedGamma4 = list(
    np.arange(99, 110 + 1, 1)
)  # SGs with fixed surrogate representation with \Gamma_4 as trivial irrep


def ir2vec(v: sp.core.add.Add, irreps: list[sp.core.symbol.Symbol]) -> np.ndarray:
    """Convert a sum of the irreps of the space group into a vector of the multiplicities
    following the order specified in the input `irreps`.

    Args:
        v (sympy.core.add.Add): sum of irreps of the space group.
        irreps (list): list of irreps indicating their order and names.

    Returns:
        numpy.ndarray: vector of multiplicities of the irrep's sum following the order
        specified in `irreps`.
    """  # noqa: E501
    y = np.array([v.coeff(irrep) for irrep in irreps])
    return y


def ebr2vec(v: sp.core.add.Add, ebrs: list[sp.core.symbol.Symbol]) -> np.ndarray:
    """Convert a sum of EBRs into a vector of multiplicities of the ebrs of the space group
    following the order stipulated in the input `ebrs`.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        ebrs (list): list of EBRs indicating their order and names.

    Returns:
        numpy.ndarray: vector of multiplicities of the ebrs following the order specified in
        the input `ebrs`.
    """  # noqa: E501
    vec = np.array([v.coeff(ebr) for ebr in ebrs])
    return vec


def ebr2irvec(
    v: sp.core.add.Add, ebrs: list[sp.core.symbol.Symbol], EBR: np.ndarray
) -> np.ndarray:
    """Convert a sum of EBRs into a vector of multiplicities of the irreps of the space group
    following the order stipulated in the input `EBR`.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        ebrs (list): list of EBRs indicating their order and names.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR.

    Returns:
        numpy.ndarray: vector of multiplicities of the irreps following the order specified
        in the input `EBR`.
    """  # noqa: E501
    vec = EBR @ np.array([v.coeff(ebr) for ebr in ebrs])
    return vec


def ebr2ir(
    v: sp.core.add.Add,
    irreps: list[sp.core.symbol.Symbol],
    ebrs: list[sp.core.symbol.Symbol],
    EBR: np.ndarray,
) -> sp.core.add.Add:
    """Convert a sum of EBRs into a sum of irreps of the space group.

    Args:
        v (sympy.core.add.Add): sum of EBRs of the space group.
        irreps (list): list of irreps indicating their order and names.
        ebrs (list): list of EBRs indicating their order and names.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR.

    Returns:
        sympy.core.add.Add: sum of irreps of the space group.
    """  # noqa: E501
    ir = irreps @ ebr2irvec(v, ebrs, EBR)
    return ir


class varArraySolutionObtainer(cp_model.CpSolverSolutionCallback):
    """A class for a solution callback used by the `CpSolver`.

    Args:
        cp_model (list): set of variables defined within the model.
    """

    def __init__(self, variables: list[cp_model.IntVar]) -> None:
        """Initiate a solution initial state.

        Args:
            variables (ortools.sat.python.cp_model.IntVar): values of the variables.
        """
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__variables = variables
        self.__solutionCount = 0
        self.__x = []

    def on_solution_callback(self) -> None:
        """Save the current solution in a list and add `1` to the solution count."""
        self.__solutionCount += 1
        self.__x.append([self.value(v) for v in self.__variables])  # type:ignore

    def solutionCount(self) -> int:
        """Returns the total number of solutions.

        Returns:
            int: Total number of solutions
        """
        return self.__solutionCount

    def solutions(self) -> np.ndarray:
        """Returns all solutions to the model in a matrix form, where each row is a particular
        solution.

        Returns:
            numpy.ndarray: Numpy array which rows represent a possible solution to the model.
        """
        return np.array(self.__x)


def searchForAllLongModes(dim: list[int], t: int) -> varArraySolutionObtainer:
    """Construct and solve a model which search for all possible set of bands with total
    dimension `t`. The output of this function can be used as a set of auxiliary modes of
    dimension `t` to input in the *enumeration algorithm*.

    Args:
        dim (list): list indicating the dimensions of the EBRs of the space group.
        t (int): positive integer indicating the total dimension of the set of bands to
        search for.

    Returns:
        varArraySolutionObtainer: returns a container which include all solutions in matrix
        form (where each row is a particular solution) and the total number of solutions.
    """

    # Initiate a model and solver.
    model = cp_model.CpModel()
    solver = cp_model.CpSolver()

    N = len(dim)

    # Creates the variables
    x = [model.new_int_var(0, 1000, f"x{i}") for i in range(N)]
    solutionObtainer = varArraySolutionObtainer(x)  # type:ignore

    # Create the constraints.
    model.add(np.dot(x, dim) == t)  # type:ignore

    # Solve.
    solver.parameters.enumerate_all_solutions = True
    solver.SearchForAllSolutions(model, solutionObtainer)

    return solutionObtainer


def vec2ebr(v: np.ndarray, ebrs: list[sp.core.add.Add]) -> list[sp.core.add.Add]:
    """Transform a vector of multiplicities of EBRs into a sum of EBRs.

    Args:
        v (numpy.ndarray): vector of multiplicities of EBRs or matrix, which each row is a
        vector of multiplicities of EBRs. Each vector must be given in the order stipulated
        by the input `ebrs`.
        ebrs (list): list of variables defining the posible EBRs of the space group and
        their order.

    Returns:
        list: list whose elements are the sum of EBRs related to the input vector or to each
        row of the input matrix.
    """
    y = v @ ebrs
    return y


def searchForAll_EBRs(
    v: np.ndarray, dim: list[int], EBR: np.ndarray, long_modes
) -> list[np.ndarray]:
    """Define and solver a model which search for all posible linear combinations of EBRs
    which, after subtracting a set of auxiliary modes until dimension `t`, have the input `v`
    as their *symmetry vector*.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps.
        dim (list): list indicating the dimensions of the EBRs of the space group.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR in the space group.
        long_modes (np.ndarray): numpy array containing the vector of multiplicities of EBRs
        to use as auxiliary modes for the model.

    Returns:
        list: list of vectors of multiplicities of EBRs which solve the problem for a certain
        set of auxiliary modes defined by the input `long_modes`.
    """
    possibleEBRs = []  # empty list that will contain all solutions

    N_ebrs = len(dim)
    N_irr = len(v)

    for i in long_modes:
        # define a model for each auxiliary mode
        model = cp_model.CpModel()
        solver = cp_model.CpSolver()
        x = [model.NewIntVar(0, 1000, f"x{i}") for i in range(N_ebrs)]
        solutionObtainer = varArraySolutionObtainer(x)  # type:ignore

        # introduce the constrains for each model using the particular auxiliary modes
        n = EBR @ i
        y = v + n
        z = y.astype(int)
        for j in range(N_irr):
            model.Add(EBR[j] @ x == z[j])

        # solve
        solver.SearchForAllSolutions(model, solutionObtainer)

        # append the solutions found for each model with the particular auxiliary modes used
        possibleEBRs.append(solutionObtainer.solutions())

    return possibleEBRs


def phys(
    x: sp.core.add.Add,
    v: np.ndarray,
    ebrs: list[sp.core.symbol.Symbol],
    EBR: np.ndarray,
    SG: int,
    irreps: list[sp.core.symbol.Symbol],
) -> int:
    """Determines if a linear combinations of EBRs given by the input `x` are transverse
    polarized to leading order in $|\mathbf{k}|$ when $\mathbf{k}, \omega \to 0$. We define
    this situation as *physical*.

    Args:
        x (sp.core.add.Add): linear combinations of EBRs of the space group.
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in the input `irreps`.
        ebrs (list): list of variables defining the EBRs of the space group in an specific
        order.
        EBR (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR.
        SG (int): number of the space group given by the notation used in the Bilbao
        Crystallographic Server.

    Returns:
        Bool: `True` if the linear combination of EBRs `x` is physical, `False` otherwise.
    """

    # differentiate between pinned and unpinned PGs and between different trivial irreps in
    # such PGs
    if SG in fixedGamma2:
        phys = (
            ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[0]) > 0
            or ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[1]) > 0
        )
    elif SG in fixedGamma4:
        phys = (
            ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[0]) > 0
            or ebr2ir(x, ebrs=ebrs, EBR=EBR, irreps=irreps).coeff(irreps[3]) > 0
        )
    else:
        phys = (
            np.sum(ebr2irvec(x, ebrs, EBR) - v - np.abs(ebr2irvec(x, ebrs, EBR) - v))
            / 2
            == 0
        )  # checks if all negative multiplicities are contained in the auxiliary modes

    return phys


def showAllResults(
    v: np.ndarray,
    t: int,
    dim: list[int],
    EBRs: np.ndarray,
    N_gamma: int,
    ebrs: list[sp.core.symbol.Symbol],
    irreps: list[sp.core.symbol.Symbol],
    SG: int,
) -> list[
    tuple[list[sp.core.add.Add], sp.core.add.Add, list[int], list[sp.core.add.Add]]
]:
    """Obtains all linear combinations of EBRs $n^{T+L}$ which can be represented by the
    *symmetry vector* `v`. In addition, it computes any necessary auxiliary modes $n^L$ to
    obtain such decompositions in EBRs, if they are physical or not, and, the surrogated
    representation at the $\Gamma$ point and zero frequency.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified in the input `irreps`.
        t (int): positive integer indicating the dimension of the auxiliary modes to search
        for.
        dim (list): list indicating the dimensions of the EBRs of the space group in the
        order specified by the input `ebrs`.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR.
        N_gamma (int): number of irreps at $\Gamma$.
        ebrs (list): list of variables defining the EBRs of the space group and their order.
        irreps (list): list of irreps of the space group indicating their order and their names.
        SG (int): number of the space group under study.

    Returns:
        list: list of tuples which first component is a list of possible lineal
        combinations of EBRs $n^{T+L}$, second component is the EBRs describing the auxiliary
        modes used $n^L$, third component is a list of bools indicating if the EBRs
        decompositions are physical (`True`) or not (`False`), and, fourth component is the
        surrogated representation at $¢$\Gamma$ and zero frequency.
    """
    # remove Gamma from the study so the algorithm is Gamma agnostic
    rv = np.delete(v, range(N_gamma))  # type: ignore
    rEBR = np.delete(EBRs, range(N_gamma), axis=0)  # type: ignore

    # compute all possible auxiliary modes with dimension t
    long_modes = searchForAllLongModes(dim, t).solutions()

    # use such auxiliary modes to search for all possible EBR decompositions of the symmetry
    # vector v
    all_EBRs = searchForAll_EBRs(rv, dim, rEBR, long_modes)

    # store the solutions found previously in a fancy way so you can trace back n^{T+L} to
    # each n^L used to compute it
    TETB_vs_LM = [
        [(all_EBRs[i] @ ebrs).tolist(), (long_modes[i] @ ebrs)]
        for i in range(len(long_modes))
        if all_EBRs[i].shape[0] > 0
    ]

    # check if the decomposition obtained is physical or not
    physical = [
        [
            phys(TETB_vs_LM[i][0][j], v, ebrs, EBRs, SG, irreps)
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    # compute the surrogate representation for each decomposition obtained
    bs = [
        [
            (
                irreps
                @ (
                    ebr2irvec(TETB_vs_LM[i][0][j], ebrs, EBRs)
                    - ebr2irvec(TETB_vs_LM[i][1], ebrs, EBRs)
                    - v
                )
            )
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    return [
        (TETB_vs_LM[i][0], TETB_vs_LM[i][1], physical[i], bs[i])
        for i in range(len(TETB_vs_LM))
    ]


def showOnlyPhysical(
    v: np.ndarray,
    t: int,
    dim: list[int],
    EBRs: np.ndarray,
    N_gamma: int,
    ebrs: list[sp.core.symbol.Symbol],
    irreps: list[sp.core.symbol.Symbol],
    SG: int,
) -> list[tuple[list[sp.core.add.Add], sp.core.add.Add, list[sp.core.add.Add]]]:
    """Obtains all *physical* linear combinations of EBRs $n^{T+L}$ which can be represented
    by the *symmetry vector* `v`. In addition, it computes any necessary auxiliary modes $n^L$
    to obtain such decompositions in EBRs, and, the surrogated representation at the $\Gamma$
    point and zero frequency. Same functionally as `showAllResults` but filtering only the
    *physical* solutions.

    Args:
        v (numpy.ndarray): vector of multiplicities of the irreps following the order
        specified by the input `irreps`.
        t (int): positive integer indicating up to what dimension of the auxiliary modes search.
        dim (list): list indicating the dimensions of the EBRs of the space group in the order
        specified by the input `ebrs`.
        EBRs (numpy.ndarray): matrix form by columns with the vectors of multiplicities of
        the irreps for each EBR.
        N_gamma (int): number of irreps at Gamma.
        ebrs (list): list of variables defining the posible EBRs of the space group and their
        order.
        irreps (list): list of irreps indicating their order and their names.
        SG (int): number of the space group under study.

    Returns:
        list: list of tuples which first component indicate the physical linear
        combination of EBRs $n^{T+L}$, second component is the EBRs representing the auxiliary
        modes used $n^L$, and, third component the surrogated representation at $\Gamma$ and
        zero frequency.
    """
    rv = np.delete(v, range(N_gamma))  # type: ignore
    rEBR = np.delete(EBRs, range(N_gamma), axis=0)  # type: ignore

    long_modes = searchForAllLongModes(dim, t).solutions()

    all_EBRs = searchForAll_EBRs(rv, dim, rEBR, long_modes)

    TETB_vs_LM = [
        [(all_EBRs[i] @ ebrs).tolist(), (long_modes[i] @ ebrs)]
        for i in range(len(long_modes))
        if all_EBRs[i].shape[0] > 0
    ]

    physical = [
        [
            phys(TETB_vs_LM[i][0][j], v, ebrs, EBRs, SG, irreps)
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    bs = [
        [
            (
                irreps
                @ (
                    ebr2irvec(TETB_vs_LM[i][0][j], ebrs, EBRs)
                    - ebr2irvec(TETB_vs_LM[i][1], ebrs, EBRs)
                    - v
                )
            )
            for j in range(len(TETB_vs_LM[i][0]))
        ]
        for i in range(len(TETB_vs_LM))
    ]

    return [
        (TETB_vs_LM[i][0][j], TETB_vs_LM[i][1], bs[i][j])
        for i in range(len(TETB_vs_LM))
        for j in range(len(TETB_vs_LM[i][0]))
        if physical[i][j] == True
    ]
