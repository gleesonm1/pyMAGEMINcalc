## pyMAGEMINcalc

Petrological calculations using the Holland (MAGEMin) thermodynamic database.

# Using MAGEMin
This package allows the user to perform calculations using the thermodynamic database of Holland et al. (2018). This is made possible through the recent development and release of MAGEMin (https://github.com/ComputationalThermodynamics; https://doi.org/10.1029/2022GC010427), with Julia functions called from Python used to run the calculations.

Before this package can be used, the Julia code used to compile and run the MAGEMin calculations, which is hosted in a separate repository, must be imported and added to Julia.
You will first need to open Python and install the PyJulia packages by running the following lines:
import julia
julia.install()

The Julia code used to perform MAGEMin calculations can then be added by:
1. Install Julia (https://julialang.org/downloads/).
2. Run Julia and add the MAGEMinCalc package via the following commands:
a. import Pkg

b. using Pkg

c. Pkg.add(url = "https://github.com/gleesonm1/MAGEMinCalc")

# Current calculations
At the moment, only a small selection of calculations are possible in pyMAGEMINcalc, but this is expected to expand rapidly in the future.

Currently, users can:
1. Search for melt liquidus temperatures.
2. Run isobaric and polybaric crystallisation models.

# Integration with PetThermoTools
pyMAGEMINcalc can be used individually or in combination with PetThermoTools (github.com/gleesonm1/PetThermoTools). While the PetThermoTools package is primarily designed to facilitate calculations using the MELTS thermodynamic models, the underlying code can also be used to call pyMAGEMINcalc calculations, allowing these calculations to be performed in parallel and significantly reducing computational time when multiple calculations are performed.

To learn how pyMAGEMINcalc can be run through PetThermoTools please see the README at github.com/gleesonm1/PetThermoTools.