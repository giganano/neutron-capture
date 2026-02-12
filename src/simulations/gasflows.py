r"""
Handles radial gas flows in these models.
"""

from .._globals import MAX_SF_RADIUS, END_TIME, M_STAR_MW
from .models.utils import get_bin_number, sinusoid
from .models.gradient import gradient
from .outflows import evoldata
from vice.toolkit.interpolation import interp_scheme_1d
from vice.milkyway.milkyway import _MAX_RADIUS_ as MAX_RADIUS # 20 kpc
from scipy.integrate import solve_ivp
import numpy as np
import warnings
import vice
from vice import ScienceWarning
import sys

class driver(interp_scheme_1d):

	def __init__(self, *args, dt = 0.01, **kwargs):
		inter_scheme_1d.__init__(self, *args, **kwargs)
		self.dt = dt
	
	def __call__(self, x):
		test = super().__call__(x)
		if test < 0:
			return 0
		else:
			return min(test, 0.01 / self.dt - 1.e-9)


class base:

	NORMTIME = 0.01 # Gyr

	def __init__(self, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		if outfilename is not None:
			self.outfile = open(outfilename, "w")
			self.outfile.write("# Time [Gyr]    ")
			self.outfile.write("Radius [kpc]    ")
			self.outfile.write("ISM radial velocity [kpc/Gyr]\n")
		else:
			self.outfile = None
		self.onset = onset
		self.dr = dr
		self.dt = dt


	def __enter__(self):
		return self


	def __exit__(self, exc_type, exc_value, exc_tb):
		self.outfile.close()
		return exc_value is None


	def write(self, time, radii, velocities):
		if self.outfile is not None:
			for i in range(len(radii)):
				self.outfile.write("%.5e\t%.5e\t%.5e\n" % (
					time, radii[i], velocities[i]))
		else: pass


	def setup(self, mw_model, **kwargs):
		vgas_alltimes = []
		n_zones = int(MAX_RADIUS / self.dr)
		times = [self.dt * i for i in range(int(END_TIME / self.dt) + 10)]
		for i in range(len(times)):
			if i > self.onset / self.dt:
				radii, vgas = self.__call__(i * self.dt, **kwargs)
				vgas_alltimes.append(vgas)
			else:
				radii = [self.dr * i for i in range(n_zones)]
				vgas = len(radii) * [0.]
				vgas_alltimes.append(vgas)
		matrix_elements_inward = []
		matrix_elements_outward = []
		for i in range(n_zones):
			areafracs_inward = []
			areafracs_outward = []
			vgas = [row[i] for row in vgas_alltimes]
			for j in range(len(times)):
				if j >= self.onset / self.dt:
					radius = i * self.dr
					if vgas[j] > 0: # outward flow
						numerator = 2 * (radius + 
							self.dr) * vgas[j] * self.NORMTIME
						numerator -= vgas[j]**2 * self.NORMTIME**2
					else: # inward flow
						numerator = vgas[j]**2 * self.NORMTIME**2
						numerator -= 2 * radius * vgas[j] * self.NORMTIME
					denominator = 2 * radius * self.dr + self.dr**2
					areafrac = numerator / denominator
					if areafrac * self.dt / self.NORMTIME > 1 - 1.e-9:
						warnings.warn("""\
Area fraction larger than 1. Consider comparing results with different \
timestep sizes to assess the impact of numerical artifacts.""", ScienceWarning)
						areafrac = 1 - 1.e-9
					elif areafrac * self.dt / self.NORMTIME < 1.e-9:
						areafrac = 1.e-9
					else: pass
					if vgas[j] > 0:
						areafracs_outward.append(areafrac)
						areafracs_inward.append(1.e-10)
					else:
						areafracs_outward.append(1.e-10)
						areafracs_inward.append(areafrac)
				else:
					areafracs_outward.append(1.e-10)
					areafracs_inward.append(1.e-10)
			matrix_elements_outward.append(
				driver(times, areafracs_outward, dt = self.dt))
			matrix_elements_inward.append(
				driver(times, areafracs_inward, dt = self.dt))
		for i in range(n_zones):
			for j in range(n_zones):
				if i - 1 == j: # inward flows
					mw_model.migration.gas[i][j] = matrix_elements_inward[i]
				elif i + 1 == j: # outward flows
					mw_model.migration.gas[i][j] = matrix_elements_outward[i]
				else:
					mw_model.migration.gas[i][j] = 0


	@staticmethod
	def area_fraction(radius, vgas, dr = 0.1):
		denominator = (radius + dr)**2 - radius**2
		if vgas > 0:
			numerator = (radius + dr)**2 - (
				radius + dr - vgas * base.NORMTIME)**2
		elif vgas < 0:
			numerator = (radius - vgas * base.NORMTIME)**2 - radius**2
		else:
			numerator = 0
		return numerator / denominator


class constant(base):

	def __init__(self, speed, onset = 1, dr = 0.1, dt = 0.01,
		outfilename = "gasvelocities.out"):
		base.__init__(self, onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.speed = speed


class constant_ifrmode(constant):

	def __init__(self, radius, *args, inward = True, **kwargs):
		self.radius = radius
		self.inward = inward
		constant.__init__(self, *args, **kwargs)

	def __call__(self, **ism_state):
		if ism_state["time"] < self.onset: return 0
		if self.speed != 0:
			if (self.inward and self.speed < 0) or (
				not self.inward and self.speed > 0):
				self.write(ism_state["time"], [self.radius], [self.speed])
				frac = self.area_fraction(self.radius, self.speed, dr = self.dr)
				if frac < 0:
					frac = 0
				elif frac > 1 - 1.0e-9:
					frac = 1 - 1.0e-9
				else: pass
				return frac
			else:
				return 0
		else:
			# without this seemingly useless if-statement, both inward and
			# outward components write to the output file, resulting in
			# duplicate entries of zero velocity.
			if self.inward: self.write(ism_state["time"], [self.radius], [0])
			return 0


class amd(base):

	def __init__(self, mw_model, beta_phi_in = 0.7, beta_phi_out = 0, onset = 1,
		dr = 0.1, dt = 0.01, outfilename = "gasvelocities.out"):
		base.__init__(self, onset = onset, dr = dr, dt = dt,
			outfilename = outfilename)
		self.mw_model = mw_model
		self.beta_phi_in = beta_phi_in
		self.beta_phi_out = beta_phi_out


class amd_ifrmode(amd):

	def __init__(self, radius, *args, inward = True, **kwargs):
		self.radius = radius
		self.inward = inward
		amd.__init__(self, *args, **kwargs)

	def __call__(self, **ism_state):
		if ism_state["time"] < self.onset: return 0
		# if callable(self.beta_phi_in):
		# 	beta_phi_in = self.beta_phi_in(self.radius, ism_state["time"])
		# else:
		# 	beta_phi_in = self.beta_phi_in
		# if callable(self.beta_phi_out):
		# 	beta_phi_out = self.beta_phi_out(self.radius, ism_state["time"])
		# else:
		# 	beta_phi_out = self.beta_phi_out
		# vgas = ism_state["ofr"] / ism_state["mgas"] * (1 - beta_phi_out)
		# vgas -= ism_state["ifr"] / ism_state["mgas"] * (1 - beta_phi_in)
		# vgas *= 1e9 # kpc/yr -> kpc/Gyr ~ km/s
		# vgas *= self.radius
		vgas = self.vgas(**ism_state)
		if vgas != 0:
			if (self.inward and vgas < 0) or (not self.inward and vgas > 0 or
				self.radius == 0):
				self.write(ism_state["time"], [self.radius], [vgas])
				frac = self.area_fraction(self.radius, vgas, dr = self.dr)
				if frac < 0:
					frac = 0
				elif frac > 1 - 1.0e-9:
					frac = 1 - 1.0e-9
				else: pass
				return frac
			else:
				return 0
		else:
			# without this seemingly useless if-statement, both inward and
			# outward components write to the output file, resulting in
			# duplicate entries of zero velocity.
			if self.inward or (not self.inward and self.radius == 0):
				self.write(ism_state["time"], [self.radius], [0])
			return 0

	def vgas(self, **ism_state):
		if ism_state["time"] < self.onset: return 0
		if callable(self.beta_phi_in):
			beta_phi_in = self.beta_phi_in(self.radius, ism_state["time"])
		else:
			beta_phi_in = self.beta_phi_in
		if callable(self.beta_phi_out):
			beta_phi_out = self.beta_phi_out(self.radius, ism_state["time"])
		else:
			beta_phi_out = self.beta_phi_out
		vgas = ism_state["ofr"] / ism_state["mgas"] * (1 - beta_phi_out)
		vgas -= ism_state["ifr"] / ism_state["mgas"] * (1 - beta_phi_in)
		vgas *= 1e9 # kpc/yr -> kpc/Gyr ~ km/s
		vgas *= self.radius
		return vgas

























































