"""Constants and unit conversions. All internal quantities are SI; names carry their unit."""

import numpy as np
import numpy.typing as npt

G0_M_S2 = 9.80665  # standard gravity, 3rd CGPM (1901)
R_UNIVERSAL_J_MOL_K = 8.314462618  # CODATA 2018
P_STANDARD_PA = 101_325.0  # standard atmosphere, USSA 1976

PA_PER_BAR = 1.0e5
PA_PER_KPA = 1.0e3
PA_PER_MPA = 1.0e6
M_PER_MM = 1.0e-3
M_PER_UM = 1.0e-6
M3_PER_L = 1.0e-3
M3_PER_CM3 = 1.0e-6
KG_PER_G = 1.0e-3


def bar_to_pa[V: (float, npt.NDArray[np.float64])](value_bar: V) -> V:
    return value_bar * PA_PER_BAR


def pa_to_bar[V: (float, npt.NDArray[np.float64])](value_pa: V) -> V:
    return value_pa / PA_PER_BAR


def mm_to_m[V: (float, npt.NDArray[np.float64])](value_mm: V) -> V:
    return value_mm * M_PER_MM


def m_to_mm[V: (float, npt.NDArray[np.float64])](value_m: V) -> V:
    return value_m / M_PER_MM
