import tempfile
import requests
import matplotlib.pyplot as pp
import brille as b, brille.plotting as bp
from pathlib import Path
from euphonic.data.interpolation import InterpolationData
from brille.euphonic import BrEu

def fetch_example_InterpolationData_object(material):
	base_url = "https://raw.githubusercontent.com/g5t/brille/master/docs"
	to_fetch = ("castep_bin",)
	tmp_dir = tempfile.TemporaryDirectory()

	for tf in to_fetch:
		r = requests.get(base_url + "/" + material + "." + tf)
		if not r.ok:
			raise Exception("Fetching {} failed with reason '{}'".format(material+"/"+tf, r.reason))
		out_path = Path(tmp_dir.name, tf)
		open(str(out_path), 'wb').write(r.content)

	idata = InterpolationData.from_castep(path=tmp_dir.name, seed=material)

	tmp_dir.cleanup()
	return idata

nacl_idata =

nacl_breu = BrEu(fetch_example_InterpolationData_object('NaCl'),
                 trellis=True, max_volume=0.001, parallel=False)
				 
