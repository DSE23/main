
import sys
sys.path.append('../')
# This makes sure the parent directory gets added to the system path

from Misc import ureg, Q_
# Imports the unit registry from the Misc folder

Extra_roll_MOM = Q_("395 ((lbf*s**2)/ft)*(ft**2)")
print(Extra_roll_MOM.to(ureg("kg*m^2")))

