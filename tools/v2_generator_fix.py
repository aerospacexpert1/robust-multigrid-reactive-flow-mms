#!/usr/bin/env python3
from pathlib import Path
p=Path('tools/generate_mms_package.py')
s=p.read_text()
s=s.replace("Ta=[1000.0,500.0,350.0][lev-1]", "Ta=[1000.0,500.0,100.0][lev-1]")
s=s.replace("invMW=YF/sp.Float(16.04)+YO/sp.Float(32.0)+YP/sp.Float(44.01)+YN/sp.Float(28.014)", "invMW=YF/sp.Float(28.01055)+YO/sp.Float(31.99880)+YP/sp.Float(44.00995)+YN/sp.Float(28.01340)")
p.write_text(s)
print('V2_GENERATOR_FIX_COMPLETE')
