# Rigorous coupled-wave analysis

There are two distinct uses of the rigorous coupled-wave analysis code: to create redistribution matrices for the matrix framework, and to define and calculate structures which are treated in their entirety by RCWA (so the matrix framework is not used). The function `RCWA` is used to create redistribution matrices, while the class `rcwa_structure` is used to define structures for RCWA calculations.

## Individual function documentation

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA(structure, size, orders, options, structpath, incidence, transmission, only_incidence_angle=False, prof_layers=None, front_or_rear='front', surf_name='', detail_layer=False, save=True, overwrite=False)`

Calculates the reflection/transmission and absorption redistribution matrices for an interface using rigorous coupled-wave analysis.

#### Parameters

* **structure** – list of Solcore Layer objects for the surface
* **size** – tuple with the vectors describing the unit cell: `((x1, y1), (x2, y2))`
* **orders** – number of RCWA orders to be used for the calculations
* **options** – user options (State object)
* **structpath** – file path where matrices will be stored or loaded from
* **incidence** – incidence medium
* **transmission** – transmission medium
* **only_incidence_angle** – if True, the calculations will only be performed for the incidence theta and phi specified in the options (rest of the matrix will be zeros). CURRENTLY DOES NOT WORK CORRECTLY
* **prof_layers** – If no profile calculations are being done, None. Otherwise a list of layers in which the profile should be calculated (front non-incidence medium layer has index 1)
* **front_or_rear** – a string, either 'front' or 'rear'; front incidence on the stack, from the incidence medium, or rear incidence on the stack, from the transmission medium.
* **surf_name** – name of the surface (to save the matrices generated).
* **detail_layer**
* **save** – whether to save the redistribution matrices (True/False), default True
* **overwrite** – whether to overwrite existing saved matrices (True/False), default False

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_structure_wl(wl, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, A_per_order, detailed_rcwa, S4_options)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_structure_wl_inkstone(wl, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, *args)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_wl(wl, geom_list, l_oc, s_oc, s_names, pol, theta, phi, widths, size, orders, phi_sym, theta_intv, phi_intv, angle_vector_0, S4_options, layer_details=False, side=1, dist=None)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_wl_inkstone(wl, geom_list, l_oc, s_oc, s_names, pol, theta, phi, widths, size, orders, phi_sym, theta_intv, phi_intv, angle_vector_0, S4_options, layer_details=False, side=1, dist=None)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_wl_prof(wl, rat_output_A, dist, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, S4_options)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.RCWA_wl_prof_inkstone(wl, rat_output_A, dist, geom_list, layers_oc, shapes_oc, s_names, pol, theta, phi, widths, size, orders, *args)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.fold_phi(phis, phi_sym)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.get_reciprocal_lattice(size, orders)`

Returns the reciprocal lattice as defined in S4 (note that this is missing a factor of 2pi compared to the standard definition).

#### Parameters

* **size** – lattice vectors in real space `((ax, ay), (bx, by))`
* **orders** – number of Fourier orders to keep

#### Returns

* reciprocal lattice (tuple)

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.initialise_S(size, orders, geom_list, mats_oc, shapes_oc, shape_mats, widths, options)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.initialise_S_inkstone(size, orders, geom_list, mats_oc, shapes_oc, shape_mats, widths)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.necessary_materials(geom_list)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_absorption_per_layer(S, n_layers, incpf)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_absorption_per_layer_order(S, n_layers, incpf)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_position_resolved(S, layer, depth, incpf)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_rt(S, n_layers)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_rt_pfbo(S, n_layers, incpf, det_l=False)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_rt_pfbo_inkstone(S, n_layers, incpf)`

---

## `class rayflare.rigorous_coupled_wave_analysis.rcwa.rcwa_structure(structure, size, options, incidence, transmission)`

Calculates the reflected, absorbed and transmitted intensity of the structure for the wavelengths and angles defined using an RCWA method implemented using the S4 package.

### Parameters

#### structure

A Solcore Structure/SolarCell object with layers and materials. Alternatively, you can supply a list which can contain any mixture of Solcore Layer objects and layers defined in one of the two following ways:

1. A list of length 4, for materials with a constant refractive index. The list entries are:

   `[width of the layer in nm, real part of refractive index (n), imaginary part of refractive index (k), geometry]`

2. A list of length 5, for materials with a wavelength-dependent refractive index. The list entries are:

   `[width of the layer in nm, wavelengths, n at these wavelengths, k at these wavelengths, geometry]`

#### size

tuple with the vectors describing the unit cell: `((x1, y1), (x2, y2))`

#### options

State object containing user options

#### incidence

semi-infinite incidence medium

#### transmission

semi-infinite transmission medium (substrate)

---

### `calculate(options)`

Calculates the reflected, absorbed and transmitted intensity of the structure for the wavelengths and angles defined.

#### Parameters

**options** – options for the calculation. The key entries are:

* wavelength: Wavelengths (in m) in which calculate the data. An array.
* theta_in: polar angle (in radians) of the incident light.
* phi_in: azimuthal angle (in radians) of the incident light.
* pol: Polarisation of the light: 's', 'p' or 'u'.
* orders: number of Fourier orders to retain in the RCWA calculation
* parallel: True of False, whether or not to run simulation on parallel
* n_jobs: if parallel, specifies how many cores are used. See joblib documentation
* A_per_order: whether or not to calculate the absorption per diffraction order
* S4_options: options passed to the S4 solver.
* RCWA_method: 'inkstone' or 'S4' (upper/lower case both accepted)

#### Returns

A dictionary with the R, A and T at the specified wavelengths and angle.

---

### `calculate_profile(options)`

It calculates the absorbed energy density within the material.

In principle this has units of [power]/[volume], but we can express it as a multiple of incoming light power density on the material, which has units [power]/[area], so that absorbed energy density has units of 1/[length].

#### Parameters

* **options** – State object containing user options

---

### `edit_geom_list(layer_index, geom_index, geom_entry)`

Allows the user to edit the geom_list of a specific layer in the structure.

#### Parameters

* **layer_index** – for which layer to edit the geom_list; 0 is the incidence medium
* **geom_index** – which entry in that layer’s geom_list to change
* **geom_entry** – new entry for the geom_list.

---

### `get_fields(layer_index, wavelength, options, extent=None, depth=1e-10, n_points=200, plot=True)`

Get the components of the E and H fields at a specific depth in a layer, over a range of x/y points. Can also plot results automatically. Uses the S4 function GetFields().

#### Parameters

* **layer_index** – index of the layer in which to get epsilon. layer 0 is the incidence medium, layer 1 is the first layer in the stack, etc.
* **wavelength** – wavelength (in nm) at which to get epsilon
* **options** – State object containing user options
* **extent** – range of x/y values in format `[[x_min, x_max], [y_min, y_max]]`. Default is None, will choose a reasonable area based on the unit cell size by default
* **depth** – depth in the layer (from the top of the layer) in nm at which to calculate the fields
* **n_points** – number of points to scan across in the x and y directions
* **plot** – plot the results (True or False, default True)

#### Returns

`xs, ys, E, H, E_mag, H_mag`

x points, y points, the complex (x, y, z) components of the E field, the complex (x, y, z) components of the H field, the magnitude of the E-field, the magnitude of the H-field. The magnitude is given by `sqrt(abs(Ex^2 + Ey^2 + Ez^2))`

---

### `get_fields_unit_cell(layer_index, wavelength, options, depth=1e-10, n_points=200)`

Get the components of the E and H fields at a specific depth in a layer, over a range of x/y points. Can also plot results automatically. Uses the S4 function GetFieldsOnGrid().

#### Parameters

* **layer_index** – index of the layer in which to get epsilon. layer 0 is the incidence medium, layer 1 is the first layer in the stack, etc.
* **wavelength** – wavelength (in nm) at which to get epsilon
* **options** – State object containing user options
* **depth** – depth in the layer (from the top of the layer) in nm at which to calculate the fields
* **n_points** – number of points to scan across in the x and y directions

#### Returns

`xs, ys, E, H, E_mag, H_mag`

x points, y points, the complex (x, y, z) components of the E field, the complex (x, y, z) components of the H field, the magnitude of the E-field, the magnitude of the H-field. The magnitude is given by `sqrt(abs(Ex^2 + Ey^2 + Ez^2))`

---

### `get_fields_z_integral(layer_index, wavelength, options, extent=None, n_points=200, plot=True)`

Get the magnitude of the E and H fields integrated over z in a layer, over a range of x/y points. Can also plot results automatically.

#### Parameters

* **layer_index** – index of the layer in which to get epsilon. layer 0 is the incidence medium, layer 1 is the first layer in the stack, etc.
* **wavelength** – wavelength (in nm) at which to get epsilon
* **options** – State object containing user options
* **extent** – range of x/y values in format `[[x_min, x_max], [y_min, y_max]]`. Default is None, will choose a reasonable area based on the unit cell size by default
* **n_points** – number of points to scan across in the x and y directions
* **plot** – plot the results (True or False, default True)

#### Returns

`xs, ys, E, H, E_mag, H_mag`

x points, y points, the (x, y, z) amplitudes squared of the E-field (|Ex|² etc.), the (x, y, z) amplitudes squared of the H-field (|Ex|² etc.), the magnitude of the E-field, the magnitude of the H-field. The magnitude is given by `sqrt(abs(Ex^2 + Ey^2 + Ez^2))`

---

### `get_fourier_epsilon(layer_index, wavelength, options, extent=None, n_points=200, plot=True)`

Get the Fourier-decomposed epsilon scanning across x-y points for some layer in the structure for the number of order specified in the options for the structure. Can also plot this automatically.

#### Parameters

* **layer_index** – index of the layer in which to get epsilon. layer 0 is the incidence medium, layer 1 is the first layer in the stack, etc.
* **wavelength** – wavelength (in nm) at which to get epsilon
* **options** – State object containing user options
* **extent** – range of x/y values in format `[[x_min, x_max], [y_min, y_max]]`. Default is None, will choose a reasonable area based on the unit cell size by default
* **n_points** – number of points to scan across in the x and y directions
* **plot** – plot the results (True or False, default True)

#### Returns

`xs, ys, a_r, a_i`

The x points, y points, and the real and imaginary parts of the dielectric function.

---

### `make_S(options, wl_ind)`

---

### `save_layer_postscript(layer_index, options, filename)`

---

### `set_size(new_size)`

Allows the user to set new basis vectors for the structure.

#### Parameters

* **new_size** – new basis vectors in the format `((x1, y1), (x2, y2))`

---

### `set_widths(new_widths)`

Allows the user to set new widths for the layers in the structure.

#### Parameters

* **new_widths** – new layer widths, in nm.

---

### `update_oc(wavelengths)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.set_incident_wave(S, s, p, options, wavelength)`

---

### `rayflare.rigorous_coupled_wave_analysis.rcwa.set_incident_wave_inkstone(S, s, p, options, wavelength)`
