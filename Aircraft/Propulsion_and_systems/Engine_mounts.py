## See tutorial/pint example with subfolders why this needs to be here!
import sys
sys.path.append('../')


from Misc import ureg, Q_

def initialise_enginemass(inp):
    global enginedrymass
    enginedrymass = inp

enginedrymass = 446
enginedrymass *= ureg.pounds

def initialise_engineIxg(inp):
    global engineIxg
    engineIxg = inp

engineIxg = Q_("84.4 inch*lbf*s**2")
print(engineIxg.to(ureg("kg*m**2")))