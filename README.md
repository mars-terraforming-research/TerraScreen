# TerraScreen
Particle Screening Tool for Mars TerraForming


The particle screening tool is a radiative transfer code suitable to simulate representative Martian atmospheric and surface temperature in the presence of natural or engineered aerosols.

# Quick start guide

## Overview and compilation

The main driver for this software package is `TerraScreen.f`. 

You may open `TerraScreen.f` in a text editor and notice various parameters,  mostly self-explanatory, that can be modified prior to compilation:
*  conrath-nu parameter  `CONRNU`  (define the vertical extend of the dust)
*  cosine of the solar zenith angle `acosz`
*  surface albedo  `ALBV`
* A set of functions (e.g. `press_to_temp()` defined at the end of the file and used to specify the initial temperature profile.
* etc ... 

The code can be run in two modes, which are selected through a flag (~Line 150): 

* When `equilibrate =.true.`, the code runs to radiative equilibrium, including a simple parameterization for convection
* When `equilibrate =.false.`, the code does a single pass through the radiation without updating the temperatures. The output can be interpreted as the radiative forcing resulting from instantenously injecting aerosols in the atmosphere.

Additionally, the code is wrapped  within a main loop used to explore sensitivity to a single paremeter. Typically, the code is ran for different visible opacities, defined within an array `tau_array`, but the code could be just as easily adapted to explore other any other parameters. 
```
    parameter (i_tau=14) 
    real*8 tau_array(i_tau)
    [...]
    do ii=1,i_tau
      TAUTOT    = tau_array(ii)
    [...]
```

After selecting `equilibrate =.true.` or  `equilibrate =.false.`The executable is built using a Makefile with: 

```bash
make clean
make
```

This will produce a `TerraScreen` executable in the current directory.  `TerraScreen` requires a single argument: the aerosols optical properties passed as an input file e.g.

`./TerraScreen data/QEXT_96IR_84_VIS_dust1.5`

Because the optical properties are parsed at runtime, the code may be executed in batch mode for different aerosols type without requiring subsequent compilation. A bash script is included as an example and may be run with: 
```
./batch_run.sh
```

## Input structure

The optical properties are passed as a file in the data directory `/data/QEXT_96IR_84_VIS_XXX` that includes 

* 3 comment lines, not parsed by the code.
* Extinction efficiency, scattering efficiency, single scattering albedo, scattering assymmetry `g` for *84 visible* bands defined by their bounds `L1` and `L2` in [um] and space-separated.
* 3 comment lines, not parsed by the code.
* Extinction efficiency, scattering efficiency, single scattering albedo, scattering assymmetry `g` for *96 infrared* bands defined by their bounds `L1` and `L2` in [um]  and space-separated.

```
[Comment line #1]
[Comment line #2]
Bin  Qext    Qscat     W0      g     L1    L2  [Comment line #3]
 1    1.449   1.355   0.935  0.495  4.50  4.36 
[...]
 84    2.193   1.346   0.614  0.909  0.25  0.24 
[Comment line #4]
IR:  [Comment line #5]
Bin  Qext    Qscat     W0      g      L1       L2 [Comment line #6]
 1    0.000   0.000   0.000  0.001  1000.00  292.68 
 [...]
 96    1.374   1.280   0.932  0.484  4.66  4.50 
```

Notes:
* the 84 visible and 96 infrared spectral bands are hard-coded in `./setspi.f90` and `./setspv.f90` respectively , not read dynamically. 
* The prefix `QEXT_96IR_84_VIS_` is expected for the input optical properties. The rest of the name `QEXT_96IR_84_VIS_`*`XXXX`* is used to name the output files.

## Output files

The output of the code are saved in the `/output/` directory as :
* `output_XXXX.txt` if the code was ran with `equilibrate =.true.`
* `static_XXXX.txt` if the code was ran with `equilibrate =.false.`

Feel free to delete the files already present in `/output` and re-generate those. 

The data structure of the `output_XXXX.txt` (ran to radiative equilibirum) and `static_XXXX.txt` (instantenous forcings) are identical:
```
Particle       , [input file]
Ncase          , [number of case ran]
dt (sols)      , [time step]
Qext 0.67um    , [reference extinction at 670 nm]
Alb sfc        , [surface albedo]
Conrath nu     , [conrath-nu parameter]
Sun Flux (W/m2), [solar flux at the top of the colum]
BWN IR (cm-1)  , [infrared bounds, in wavenumber]
BWN VIS (cm-1) , [visible bounds , in wavenumber]
================
wavelenght [um], [Center wavelenght for each band VIS>IR]  
Qext           , [Exctinction efficiency for each band VIS>IR]  
Qscat          , [Scattering efficiency for each band VIS>IR] 
G factor       , [Assymetry factor for each band VIS>IR]
================
P [mbar]       ,  [layer mid-ppoints pressures in mbar]
Temp initial[K],  [initial temperature [profile]
================

it   ,  tau  , Tsfc , alb , OLR  , ASR ,  NET top , NET bot , T lev[k], OLR [nIR],ASR [nVIS],OLR [nIR],ASR [nVIS],SFC DIR[nIR],SFC DVIS[nVIS]
[...]
```
With
* `it`: number of iterations to convergence, if applicable
* `tau`: visible opacity at 670nm
* `Tsfc`: surface temperature [K]
* `alb` : top of the atmosphere (TOA) albedo = upward vis/downward vis 
* `OLR` : outgoing longwave radiation, integrated value [W/m2]
* `ASR` : absorbed shortwave radiation (down VIS- up VIS), integrated value [W/m2]
* `NET top`: radiative budget at the TOA: down IR +down VIS - up IR -up VIS [W/m2]
* `NET bot`: radiative budget at the bottom of atmosphere (BOA): down IR +down VIS - up IR -up VIS [W/m2]
* `T lev[k]`: temperature at the model's pressure levels
* `OLR [nIR]`: Outgoing longwave radiation in the 96 IR bands [W/m2/um]
* `ASR [nVIS]`: Absorbed solar radiation (down VIS- up VIS) at the TOA in the 84 bands [W/m2/um]
* `SFC DIR[nIR]`: Downward radiation at the BOA for the 96 IR bands [W/m2/um]
* `SFC DVIS[nVIS]`: Downward radiation at the BOA for the 84 VIS bands [W/m2/um]

Notes:
* The widths for the IR (resp. VIS) bands can be calculated as `DWNI(N)= BWNI(N+1)-BWNI(N)` (resp. `DWNV(N)= BWNV(N+1)-BWNV(N)`) in [cm-1]
* The center wavenumbers for the IR (resp. VIS) bands can be calculated as  `WNI(N)=BWNI(N)+DWNI(N)/2` (resp. `WNV(N)=BWNV(N)+DWNV(N)/2`) in [cm-1]
* The conversion from wavenumber to wavelenght for either the centers or the edges of the IR (resp. VIS) bands are `WLI=10**4/WNI` (resp. `WLV=10**4/WNV`) in [um]
* `wavelenght`, `Qext`, `Qscat` and `G` in the table above are provided for both the VIS and IR bands by increasing *wavelength*. Therefore, those arrays would need to be inverted from *right* to *left* to match the `BWN IR` and `BWN VIS` bands.
## Requirements
This code  requires a Fortran compiler (e.g. gfortran). Some optional pre-processing and analysis scripts are provided in Python.

## Credits
This project is based on the [NASA Ames Legacy GCM radiation code](https://github.com/nasa/legacy-mars-global-climate-model). This is an initial commit. Further update of the branch will be provided as patch to the NASA Ames Legacy GCM radiation code to properly credit the authors.  

## License
-------
This project is licensed under the MIT License - see the LICENSE file for details.

