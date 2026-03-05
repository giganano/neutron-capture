
import numpy as np
import vice
import os

def subsample_stars(stars, n = 1000):
	r"""
	Subsample from the stellar populations in the model based on their
	individual masses. Individually sampled points can be interpreted as single
	stars as opposed to single stellar populations.
	"""
	assert isinstance(stars, vice.dataframe)
	assert "mass" in stars.keys()
	assert isinstance(n, int)
	total_mass = sum(stars["mass"])
	mass_fracs = [m / total_mass for m in stars["mass"]]
	indices = np.random.choice(len(stars["mass"]), size = n, p = mass_fracs)
	sub = {}
	for key in stars.keys():
		sub[key] = [stars[key][i] for i in indices]
	return vice.dataframe(sub)


def get_velocity_profile(output, lookback):
	raw = np.genfromtxt("%s_gasvelocities.out" % (output.name))
	time = output.zones["zone0"].history["time"][-1] - lookback
	diff = [abs(_ - time) for _ in output.zones["zone0"].history["time"]]
	idx = diff.index(min(diff))
	time = output.zones["zone0"].history["time"][idx]
	radii = []
	vgas = []
	for i in range(len(raw)):
		if raw[i][0] == time:
			radii.append(raw[i][1])
			vgas.append(raw[i][2])
		else: pass
	return [radii, vgas]


def get_velocity_evolution(output, radius, zone_width = 0.1):
	if os.path.exists("%s_gasvelocities.out" % (output.name)):
		raw = np.genfromtxt("%s_gasvelocities.out" % (output.name))
	else:
		lookback = output.zones["zone0"].history["lookback"]
		vgas = len(lookback) * [0.]
		return [lookback, vgas]
	zone = int(radius / zone_width)
	radius = zone * zone_width # use inner edge for sake of lookup in file
	time = []
	vgas = []
	for i in range(len(raw)):
		if raw[i][1] == radius:
			time.append(raw[i][0])
			vgas.append(raw[i][2])
		else: pass
	lookback = [time[-1] - t for t in time]
	return [lookback, vgas]


