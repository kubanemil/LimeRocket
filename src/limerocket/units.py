"""Constants and unit conversions. All internal quantities are SI; names carry their unit."""

G0_M_S2 = 9.80665  # standard gravity, 3rd CGPM (1901)
R_UNIVERSAL_J_MOL_K = 8.314462618  # CODATA 2018
P_STANDARD_PA = 101_325.0  # standard atmosphere, USSA 1976

PA_PER_BAR = 1.0e5
PA_PER_PSI = 6_894.757293168
M_PER_MM = 1.0e-3
M_PER_FT = 0.3048
M3_PER_L = 1.0e-3


def bar_to_pa(value_bar: float) -> float:
    return value_bar * PA_PER_BAR


def pa_to_bar(value_pa: float) -> float:
    return value_pa / PA_PER_BAR


def psi_to_pa(value_psi: float) -> float:
    return value_psi * PA_PER_PSI


def mm_to_m(value_mm: float) -> float:
    return value_mm * M_PER_MM


def m_to_mm(value_m: float) -> float:
    return value_m / M_PER_MM
