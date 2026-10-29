"""Export all independent numerical fixtures and a prespecified main subset for R."""
from prepare import V,dump
import numpy as np
paths=[V/'results/tests/numerical_fixture.npz']
paths+=sorted((V/'results/fixtures').glob('*r0*'))
out=V/'work/r_fixtures';out.mkdir(parents=True,exist_ok=True)
for path in paths:
    with np.load(path) as d:dump(out/(path.stem+'.json'),{k:d[k].tolist() for k in d.files})
print(len(paths),'R fixtures exported')
