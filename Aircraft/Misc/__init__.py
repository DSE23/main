"""
IMPORTANT FILE
DO NOT EDIT/DELETE
SEE THE FOLLOWING LINK WHY THIS NEEDS TO BE HERE:
    https://pint.readthedocs.io/en/latest/tutorial.html
"""
from pint import UnitRegistry, set_application_registry
ureg = UnitRegistry()
set_application_registry(ureg)
Q_ = ureg.Quantity
