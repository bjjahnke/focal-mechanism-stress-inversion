"""Infer stress from earthquake focal mechanisms when the true fault plane is unknown.

Pipeline:
    1. compute_nodal_planes   - add the second nodal plane and P/T/N axes
    2. run_stress_inversions  - repeated iterative stress inversions (Vavrycuk, 2014)
    3. summarize_inversions   - group runs into stress solutions and validate them
"""

__version__ = "1.0.0"

from .catalog import CatalogError, load_catalog, validate_catalog  # noqa: E402
from .iterative_inversion import run_stress_inversions  # noqa: E402
from .linear_inversion import fault_instability, invert_stress  # noqa: E402
from .nodal_planes import compute_nodal_planes  # noqa: E402
from .summarize import summarize_inversions  # noqa: E402

__all__ = [
    "CatalogError", "load_catalog", "validate_catalog", "compute_nodal_planes",
    "invert_stress", "fault_instability", "run_stress_inversions", "summarize_inversions",
]
