# <codecell> Import Required Modules and env vars
import os, sys

SCRIPTS = os.environ['SCRIPTS']
sys.path.append(SCRIPTS)

import myUQlib   # Import myUQlib
cwd = os.getcwd() + "/"

# <codecell> Polynomail order, number of dimensions, hyperbolic trunc factor
n = 3
d = 20
q = 0.1

# <codecell> calling uqInfo function
myUQlib.uqInfo(n, d, q)
