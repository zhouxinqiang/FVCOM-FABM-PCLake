## The **FABM-PCLake** model 

FABM-PCLake is an open source aquatic ecosystem model originally developed by Fenjuan Hu et al., at Aarhus University, Denmark, and is a further development of the original PCLake model by Jan Janse (1997). The model describes interactions between multiple trophic levels, including piscivorous, zooplanktivorous and benthivorous fish, zooplankton, zoobenthos, phytoplankton and rooted macrophytes. The model accounts for oxygen dynamics and a fully closed nutrient cycle for nitrogen and phosphorus. In contrast to the original PCLake model, FABM-PCLake enable coupling to hydrodynamic models (e.g., GOTM and GETM), and a range of additional features, including bottom-shear-dependent resuspension, a choice between multiple light-limitation functions for primary producers, a flexible configuration of the conceptual model, and coupling to other biogeochemical models available through FABM.

FABM-PCLake was originally published in: 

**FABM-PCLake – linking aquatic ecology with hydrodynamics**, Geoscientific Model Development 9: 2271-2278,
by Hu, F., Bolding, K., Bruggeman, J., Jeppesen, E., Flindt, M. R., van Gerven, L., Janse, J. H., Janssen, A. B. G., Kuiper, J. J., Mooij, W. M., and Trolle, D. 2016. 

The original PCLake model was published in:
 
**A model of nutrient dynamics in shallow lakes in relation to multiple stable states**, Hydrobiologia 342/343: 1-8, by Janse, J. 1997.



###### If you have questions and suggestions regarding this model, please contact the AU developer team:
* Fenjuan Hu: fenjuan@bios.au.dk 
* Dennis Trolle: trolle@bios.au.dk                                            
* Karsten Bolding: bolding@bios.au.dk
* Anders Nielsen: an@bios.au.dk



## Revision history
Original implementation of FABM-PCLake: xx Dec. 2016, by Fenjuan Hu and AU developer team


## Module overview
FABM-PCLake is comprised of 11 Fortran modules, including:

* model_library.f90
* auxiliary.f90
* utility.f90
* abiotic_water.f90
* abiotic_sediment.f90
* phytoplankton_water.f90
* phytoplankton_sediment.f90
* macrophytes.f90
* zooplankton.f90
* zoobenthos.f90
* fish.f90

The overall tasks of individual modules are described in brief below. 
Each module is described in greater detail on the Wiki pages (link).

### model_library
This module is a FABM wrapper module for initialising all the modules of FABM-PCLake.

### auxiliary

This modules handles general resuspension and sedimentation processes, which are used by several other modules. The module also include sediment burial processes. There are no local state variables registered in this module.


### utility

This module includes 3 types of temperature functions that are used throughout the FABM-PCLake model: 

* uFunTmAbio for all abiotic processes 
* uFunTmBio for biological processes 
* uFunTmVeg for macrophytes processes


### abiotic_water

This module describes the state variables which are related to abiotic processes and state variables in the water column, including: inorganic matter (IM), organic matter (POM and DOM), dissolved nutrients (ammonium, nitrate, phosphate and dissolved silica dioxide), immobilized phosphorus (absorbed phosphorus) and dissolved oxygen. 

The processes described in this module include mineralisation, nitrification, denitrification, phosphorus sorption and oxygen reaeration.  


### abiotic_sediment

This module describes the processes and state variables which are related to abiotic processes in the sediment, including: inorganic matter (IM), organic matter (POM and DOM), dissolved nutrients (ammonia, nitrate and phosphate) immobilized phosphorus (adsorbed phosphorus).

The processes described in this module include diffusion across sediment-water interface, mineralisation, nitrification, denitrification, phosphorus sorption. The module also describes the sediment oxic layer fraction, which is an important diagnostic variable used in several process-descriptions, e.g. relating to sorption of phosphorus. Settling and resuspension processes are described in the auxiliary module.

### phytoplankton_water

This module describes the processes and state variables related to phytoplankton in the water column, including: diatoms, green algae and cyanobacteria (blue-green algae). Each phytoplankton group is accounted for in three elements including dry-weight, nitrogen and phosphorus. The silica concentration in diatoms is not a state variables, but a diagnostic variable, as the model assumes a fixed Si/DW ratio for diatoms (i.e. 0.1).

The processes described in this module include assimilation (primary production), nutrient uptake, respiration,
excretion and mortality. The zooplankton module influence phytoplankton through grazing. 

### phytoplankton_sediment

This module describes the processes and state variables related to phytoplankton which have settled to the bottom sediments. Three groups of phytoplankton are described, including diatoms, green algae and cyanobacteria (blue-green algae). 

The processes described in this module include respiration, excretion and mortality. Settling and resuspension processes are described in the auxiliary module.

### macrophytes

This module describes the processes and state variables related to submerged macrophytes, and is implemented as a benthic module. Macrophytes are accounted for in three elements including dry-weight, nitrogen and phosphorus, and is further subdivided into a (dynamic) shoot and root fraction.

The processes described in this module include assimilation (primary production), nutrient updake, respiration, excretion, mortality and migration. The latter relates to the user option of including a colonisation rate for macrophytes, e.g., originating from the surrounding catchment. The macrophyte module influence resuspension described in the auxiliary module. The macrophyte module also contain the macrophytes coverage percentage, which is a key diagnostic variable used by the fish module.

### zooplankton

This module describes the processes and state variables relating to zooplankton. Zooplankton are accounted for in three elements including dry-weight, nitrogen and phosphorus. 

The processes described in this module include phytoplankton grazing and assimilation, respiration, excretion and mortality. The fish module influence zooplankton through predation. 

### zoobenthos

This module describes the processes and state variables relating to zoobenthos. Zoobenthos are accounted for in three elements including dry-weight, nitrogen and phosphorus. Zoobenthos feed on detritus in the sediments, and on settled phytoplankton (in the phytoplankton_sediment module). 

The processes described in this module include consumption and assimilation, migration, respiration, excretion and mortality. The fish module influence zoobenthos through predation. 

### fish

This module describes the processes and state variables relating to fish, including zooplanktivorous  (sDFiJv, sNFiJv, sPFiAd), benthivorous (sDFiAd, sNFiAd, sPFiAd) and piscivorous fish (sDPisc).

The processes described in this module include migration, reproduction (as part of adult benthivorous fish become juvenile zooplanktivorous fish), aging (as part of juvenile zooplanktivorous fish become adult benthivorous fish), assimilation (predation on zooplankton, zoobenthos or fish), respiration, excretion and mortality. In addition, the model user can choose to include a "harvest" or "stocking" rate for any of the individual fish groups.