 # FABM-MEDUSA

FABM-MEDUSA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. FABM-MEDUSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License provided in COPYING for more details.

Copyright 2020 Plymouth Marine Laboratory.

 ## About FABM-MEDUSA
FABM-MEDUSA is a modular implementation of the Model of Ecosystem Dynamics, nutrient Utilization, Sequestration and Acidification (MEDUSA) within FABM (Framework for Aquatic Biogeochemical Models; http://fabm.net) (Bruggeman and Bolding, 2014).

This implementation is based on the model description in Yool, A., Popova, E. E., and Anderson, T. R.: MEDUSA-2.0: an intermediate complexity biogeochemical model of the marine carbon cycle for climate change and ocean acidification studies, Geosci. Model Dev., 6, 1767-1811, https://doi.org/10.5194/gmd-6-1767-2013, 2013, and reproduces model functionality corresponding to NEMO revision 10020, UKMO branch /branches/UKMO/dev_r5518_GO6_package, with key_roam defined.

 ## Testcases
An example configuration file fabm-medusa-original.yaml with a default parameter set and rather arbitrary initial conditions is supplied. For simulations driven by GOTM (https://gotm.net/), sample fabm_input.nml and output.yaml files are provided. For running the model with NEMO (https://www.nemo-ocean.eu/), an output control file iodef.xml is included.

 ## Implementation notes
Original comments from NEMO-MEDUSA code are included in this version of the model, mostly in adapted form, to facilitate transition between the code bases, augmented with FABM-specific information as necessary.

**Fast-sinking detritus.** The remineralisation method currently implemented in FABM is the Ballast model (corresponding to jexport = 1, iball = 1 in MEDUSA), which couples C, N and Fe remineralisation to that of particulate Si and CaCO3. The ballast-sans-ballast model (jexport = 1, iball /= 1), Martin et al (1987) and Henson et al (2011) (jexport = 2 and 3, respectively) models are not yet implemented in FABM-MEDUSA.

**Carbonate system.** PML carbonate system is re-implemented in FABM-MEDUSA (see description in Butenschoen et al, Geosci. Model Dev., 9, 1293–1339, https://doi.org/10.5194/gmd-9-1293-2016, and references therein). MOCSY carbonate system routines (activated with key_mocsy in NEMO-MEDUSA) are not yet implemented in FABM-MEDUSA.

**Organic matter treatment at the seafloor.** In contrast to combinations of jfdfate (fate of fast detritus at the seafloor), jorgben and jinorgben (fate of organic material at the seafloor), FABM-MEDUSA implements a "seafloor" switch:

   For fast detritus (in medusa_fast_detritus):
  - `seafloor.eq.1`: (`jfdfate.eq.0 .and. jorgben.eq.0`) - immediate remineralisation (no benthos)
  - `seafloor.eq.2`: (`jfdfate.eq.1 .and. jorgben.eq.0`) - conversion of fast C into slow C (no benthos)
  - `seafloor.eq.3`: (`jfdfate.eq.0 .and. jorgben.eq.1`) - conversion of fast detritus into benthic pool

   For slow detritus (in medusa_pelagic):
  - `seafloor.eq.3`: (`jorgben.eq.1`) - benthic module is registered and the slow detritus enters the benthic pool
  - `seafloor.ne.3`: (`jorgben.eq.0`) - benthic module not registered, slow detritus stays in pelagic.

 ## Contact

For any requests regarding FABM-MEDUSA contact Gennadi Lessin (gle@pml.ac.uk).

For questions related to FABM contact Jorn Bruggeman (jorn@bolding-bruggeman.com) and visit https://fabm.net.

For questions on MEDUSA process descriptions contact Andrew Yool (axy@noc.ac.uk).
  
