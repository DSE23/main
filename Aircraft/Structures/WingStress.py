"""                  
Name: Wing_Stress_Calculations 
Department: Structures 
Last updated: 08/06/2018 12:38 by Boris 
"""

from Structures import Inertia


def Normal_stress_due_to_bending(zs, cs, y): # Normal stress due to bending
    denominator_inertia_term = Inertia.Ixx*Inertia.Iyy-Inertia.Ixy**2
    inertia_term_1 = (Inertia.Iyy*y-Inertia.Ixy*cs)/denominator_inertia_term
    inertia_term_2 = (Inertia.Ixx*cs-Inertia.Ixy*y)/denominator_inertia_term
    sigma_zs = Mx*inertia_term_1 + My*inertia_term_2
    return sigma_zs #Gives the normal stress function for a given span zs, and x- and y- coordinate


