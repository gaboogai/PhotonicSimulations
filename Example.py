import numpy as np
from solcore import si, material
from solcore.structure import Layer
from solcore.solar_cell import SolarCell

from rayflare.rigorous_coupled_wave_analysis import rcwa_structure
from rayflare.options import default_options

import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler

import warnings
warnings.filterwarnings('ignore')

# ---------- SECTION 2 ------------------------------------------------

InAlP_hole_barrier = material('AlInP')(Al=0.5)
GaAs_pn_junction = material('GaAs')()
InGaP_e_barrier = material('GaInP')(In=0.5)
Ag = material('Ag')()
Au = material('Au')()
SiN = material('Si3N4')()
Air = material('Air')()

wavelengths = np.linspace(303, 1000, 50) * 1e-9

x = 500

size = ((x, 0), (x / 2, np.sin(np.pi / 3) * x))

options = default_options()
options.wavelength = wavelengths
options.orders = 60
options.RCWA_method = 'S4'

# ---------- SECTION 3 ------------------------------------------------

InAlP_hole_barrier = material('AlInP')(Al=0.5)
GaAs_pn_junction = material('GaAs')()
InGaP_e_barrier = material('GaInP')(In=0.5)
Ag = material('Ag')()
Au = material('Au')()
SiN = material('Si3N4')()
Air = material('Air')()

wavelengths = np.linspace(303, 1000, 50) * 1e-9 #wavelengths, actually set below

x = 500 

size = ((x, 0), (x / 2, np.sin(np.pi / 3) * x))

options = default_options()
options.wavelength = wavelengths
options.orders = 60 #number of fourier orders? not quite sure what this means 
options.RCWA_method = 'S4'

#--------- Section 3 --------------------------------------------------

grating_circles =[{'type': 'circle', 'mat': Ag, 'center': (0, 0), 'radius': 115}]

grating_squares = [{'type': 'rectangle', 'mat': Ag, 'center': (0, 0),
                    'halfwidths': [115, 115], 'angle': 20}]

grating_ellipse = [{'type': 'ellipse', 'mat': Ag, 'center': (0, 0),
                    'halfwidths': [150, 110], 'angle': 20}]

grating_polygon = [{'type': 'polygon', 'mat': Ag, 'center': (0, 0), 'angle': 0,
                   'vertices': ((150, 0), (0, 150), (-150, 0))}]

grating_circle_polygon = [{'type': 'circle', 'mat': Ag, 'center': (0, 0), 'radius': 115},
                         {'type': 'polygon', 'mat': Au, 'center': (200, 150), 'angle': -20,
                   'vertices': ((150, 0), (0, 150), (-150, 0))}]

grating_list = [grating_circles, grating_squares, grating_ellipse,
                grating_polygon, grating_circle_polygon, None]

#--------- Section 4 --------------------------------------------------

GaAs_absorption = np.zeros((len(grating_list), len(wavelengths)))
grating_absorption = np.zeros((len(grating_list), len(wavelengths)))

for i1, geometry in enumerate(grating_list):
    solar_cell = SolarCell([Layer(material=InGaP_e_barrier, width=si('19nm')),
                                Layer(material=GaAs_pn_junction, width=si('85nm')),
                                Layer(material=InAlP_hole_barrier, width=si('19nm')),
                            Layer(material=SiN, width=si('100nm'), geometry=geometry)])

    S4_setup = rcwa_structure(solar_cell, size=size, options=options,
                              incidence=Air, transmission=Ag)
    RAT = S4_setup.calculate(options)
    GaAs_absorption[i1] = RAT['A_per_layer'][:, 1]
    grating_absorption[i1] = RAT['A_per_layer'][:, 3]
    pattern = S4_setup.get_fourier_epsilon(layer_index=4, wavelength=500, options=options)

#--------- Section 4 --------------------------------------------------

pal = sns.color_palette('husl', len(grating_list))

cols = cycler('color', pal)
plt.rcParams.update({'axes.prop_cycle': cols})

plt.figure()
plt.plot(wavelengths*1e9, GaAs_absorption.T)
plt.legend(['Circles', 'Squares', 'Ellipse',
            'Polygon (triangle)', 'Circles + triangles', 'Planar SiN'])
plt.plot(wavelengths*1e9, grating_absorption.T, '--')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Absorption in 80 nm GaAs/grating')
plt.autoscale(tight=True)
plt.show()