from pint import UnitRegistry

unit = UnitRegistry()

def initialise_enginemass(inp):
    global enginedrymass
    enginedrymass = inp

enginedrymass = 446
enginedrymass *= unit.pounds

def initialise_engineIxg(inp):
    global engineIxg
    engineIxg = inp

engineIxg = unit("84.4 inch*lb*s**2")