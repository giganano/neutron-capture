r"""
The diskmodel objects employed in the Johnson et al. (2021) study.
"""

try:
	ModuleNotFoundError
except NameError:
	ModuleNotFoundError = ImportError
try:
	import vice
except (ModuleNotFoundError, ImportError):
	raise ModuleNotFoundError("Could not import VICE.")
if vice.version[:2] < (1, 2):
	raise RuntimeError("""VICE version >= 1.2.0 is required to produce \
Johnson et al. (2021) figures. Current: %s""" % (vice.__version__))
else: pass
from vice.yields.presets import JW20
from vice.toolkit import hydrodisk
vice.yields.sneia.settings['fe'] *= 10**0.1
from .._globals import END_TIME, MAX_SF_RADIUS, ZONE_WIDTH
from . import gasflows
from . import migration
from . import models
from .models.utils import get_bin_number, interpolate
from .models.gradient import gradient
from .models import mergers
from . import inputs
from . import sfe
import math as m
import sys


class diskmodel(vice.milkyway):

	r"""
	A milkyway object tuned to the Johnson et al. (2021) models specifically.

	Parameters
	----------
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.
	name : ``str`` [default : "diskmodel"]
		The name of the model; the output will be stored in a directory under
		this name with a ".vice" extension.
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the star formation history.
		Allowed values:

		- "static"
		- "insideout"
		- "lateburst"
		- "outerburst"

	verbose : ``bool`` [default : True]
		Whether or not the run the models with verbose output.
	migration_mode : ``str`` [default : "diffusion"]
		A keyword denoting the time-dependence of stellar migration.
		Allowed values:

		- "diffusion"
		- "linear"
		- "sudden"
		- "post-process"

	kwargs : varying types
		Other keyword arguments to pass ``vice.milkyway``.

	Attributes and functionality are inherited from ``vice.milkyway``.
	"""

	def __init__(self, zone_width = 0.1, name = "diskmodel", spec = "static",
		verbose = True, migration_mode = "diffusion", **kwargs):
		super().__init__(zone_width = zone_width, name = name,
			verbose = verbose, **kwargs)
		if self.zone_width <= 0.2 and self.dt <= 0.02 and self.n_stars >= 6:
			Nstars = 3102519
		else:
			Nstars = 2 * int(MAX_SF_RADIUS / zone_width * END_TIME / self.dt *
				self.n_stars)
		self.migration.stars = migration.diskmigration(self.annuli,
			N = Nstars, mode = migration_mode,
			filename = "%s_analogdata.out" % (name))
		kwargs = {
			"spec": spec,
			"zone_width": zone_width
		}
		if spec in ["expifr", "expifr_gse"]:
			self.mode = "ifr"
			self.evolution = accretion_history(**kwargs)
		else:
			self.mode = "sfr"
			self.evolution = star_formation_history(**kwargs)
		self.dt = timestep_size
		self.elements = elements

		for i in range(self.n_zones):
			area = m.pi * ZONE_WIDTH**2 * ((i + 1)**2 - i**2)
			self.zones[i].Mg0 = 0
			self.zones[i].tau_star = sfe.sfe(area, mode = self.mode)
			self.zones[i].RIa = "exp"
			self.zones[i].tau_ia = 1.5

		for i in range(self.n_zones):
			self.zones[i].Zin = {}
			for elem in self.zones[i].elements:
				if spec == "expifr_gse":
					self.zones[i].Zin[elem] = mergers.Zin_with_GSE(
						zone_width * (i + 0.5),
						elem,
						dr = zone_width,
						dt = self.dt)
				else:
					self.zones[i].Zin[elem] = mergers.Zin_CGM(elem)

		if inputs.RADIAL_GAS_FLOWS is not None:
			kwargs = {
				"onset": inputs.RADIAL_GAS_FLOW_ONSET,
				"dr": zone_width,
				"dt": self.dt
			}
			if self.mode == "ifr":
				self.migration.gas.callback = True
				kwargs["outfilename"] = None
			else:
				kwargs["outfilename"] = "%s_gasvelocities.out" % (self.name)
			callkwargs = {}
			if inputs.RADIAL_GAS_FLOWS == "constant":
				if self.mode == "ifr":
					for i in range(self.n_zones):
						for j in range(self.n_zones):
							if spec == "expifr":
								obj = gasflows.constant_ifrmode
							elif spec == "expifr_gse":
								raise ValueError("Need to set up gas flows.")
							else:
								raise ValueError("Bruh.")
							if abs(i - j) == 1:
								self.migration.stars[i][j] = obj(
									i * zone_width,
									inputs.RADIAL_GAS_FLOW_SPEED,
									inward = i < j,
									**kwargs)
							else: pass
				else:
					self.radialflow = gasflows.constant(
						inputs.RADIAL_GAS_FLOW_SPEED, **kwargs)
			else:
				raise ValueError(
					"Unrecognized radial gas flow setting: %s" % (
						inputs.RADIAL_GAS_FLOWS))

			if self.mode == "ifr":
				self.outfile = open("%s_gasvelocities.out" % (self.name), 'w')
				self.outfile.write("# Time [Gyr]    ")
				self.outfile.write("Radius [kpc]    ")
				self.outfile.write("ISM radial velocity [kpc/Gyr]\n")
				for i in range(self.n_zones):
					for j in range(self.n_zones):
						if isinstance(self.migration.gas[i][j], gasflows.base):
							self.migration.gas[i][j].outfile = self.outfile
						else: pass
			else:
				self.radialflow.setup(self, **callkwargs)
		else: pass


	def run(self, *args, **kwargs):
		out = super().run(*args, **kwargs)
		self.migration.stars.close_file()
		if self.mode == "ifr" and inputs.RADIAL_GAS_FLOWS is not None:
			self.outfile.close()
		else: pass
		return out

	@classmethod
	def from_config(cls, config, **kwargs):
		r"""
		Obtain a ``diskmodel`` object with the parameters encoded into a
		``config`` object.

		**Signature**: diskmodel.from_config(config, **kwargs)

		Parameters
		----------
		config : ``config``
			The ``config`` object with the parameters encoded as attributes.
			See src/simulations/config.py.
		**kwargs : varying types
			Additional keyword arguments to pass to ``diskmodel.__init__``.

		Returns
		-------
		model : ``diskmodel``
			The ``diskmodel`` object with the proper settings.
		"""
		model = cls(zone_width = config.zone_width,
			timestep_size = config.timestep_size, elements = config.elements,
			**kwargs)
		model.n_stars = config.star_particle_density
		model.bins = config.bins
		return model


class evol_spec:

	def __init__(self, spec = "static", zone_width = 0.1):
		self._radii = []
		self._evol = []
		i = 0
		max_radius = 20 # kpc, defined by ``vice.milkyway`` object.
		while (i + 1) * zone_width < max_radius:
			self._radii.append((i + 0.5) * zone_width)
			self._evol.append({
					"oscil":		models.insideout_oscil,
					"static": 		models.static,
					"insideout": 	models.insideout,
					"lateburst": 	models.lateburst,
					"outerburst": 	models.outerburst,
					"expifr": 		models.expifr,
					"expifr_gse": 	models.expifr_with_GSE
				}[spec.lower()]((i + 0.5) * zone_width, dr = zone_width))
			i += 1


class star_formation_history(evol_spec):

	r"""
	The star formation history (SFH) of the model galaxy. This object will be
	used as the ``evolution`` attribute of the ``diskmodel``.

	Parameters
	----------
	spec : ``str`` [default : "static"]
		A keyword denoting the time-dependence of the SFH.
	zone_width : ``float`` [default : 0.1]
		The width of each annulus in kpc.

	Calling
	-------
	- Parameters

		radius : ``float``
			Galactocentric radius in kpc.
		time : ``float``
			Simulation time in Gyr.
	"""

	def __call__(self, radius, time):
		# The milkyway object will always call this with a radius in the
		# self._radii array, but this ensures a continuous function of radius
		if radius > MAX_SF_RADIUS:
			return 0
		else:
			idx = get_bin_number(self._radii, radius)
			if idx != -1:
				result = gradient(radius) * interpolate(self._radii[idx],
					self._evol[idx](time), self._radii[idx + 1],
					self._evol[idx + 1](time), radius)
			else:
				result = gradient(radius) * interpolate(self._radii[-2],
					self._evol[-2](time), self._radii[-1], self._evol[-1](time),
					radius)
			if result < 0: result = 0
			return result


class accretion_history(evol_spec):

	def __call__(self, radius, time):
		if radius > MAX_SF_RADIUS:
			return 0
		else:
			idx = get_bin_number(self._radii, radius)
			if idx != -1:
				result = interpolate(self._radii[idx],
					self._evol[idx](time), self._radii[idx + 1],
					self._evol[idx + 1](time), radius)
			else:
				result = interpolate(self._radii[-2],
					self._evol[-2](time), self._radii[-1], self._evol[-1](time),
					radius)
			if result < 0: result = 0
			return result

