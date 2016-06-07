from __future__ import division,absolute_import,unicode_literals,print_function
from six import PY2
if PY2:
    from pathlib2 import Path
else:
    from pathlib import Path
#%% apply future to each module
flist = Path(__file__).parent.glob('*.py')
flist = [f for f in flist if f!='__init__.py']
__all__=flist

for f in flist:
    from . import f #thereby applying __future__
