def prof(fn):
    return fn

import __builtin__
if not "profile" in __builtin__.__dict__:
    __builtin__.__dict__['profile'] = prof

from common import *
from mesh import *
from shapeset_h1 import *
from space_h1 import *
from discrete import *
from function import *
from view import *
from integrals_h1 import *
from solvers import *
