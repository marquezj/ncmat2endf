# ncmat2endf
Python script to generate TSL ENDF files from NCrystal

## Description

This is a script for creating a set of ENDF-6 thermal scattering files from a .ncmat
file. Parameters for the ENDF-6 file can be defined with command line arguments or
changing the endf_defaults dictionary in the script. 

The script allows to handle multiple temperatures in one ENDF-6 file, but this is not
recommended, because NCrystal computes an optimal (alpha, beta) grid for each material
and temperature.

Ths script uses the endf-parserpy package from IAEA to format and check the syntaxis of
the ENDF-6 file.

## Installation 

This script requires the [`endf-parserpy`](https://endf-parserpy.readthedocs.io/en/latest/) package developed by IAEA. This
package can be installed using `pip`:

```
pip install endf-parserpy --upgrade
```

## Usage

The script can be used from the command line:

```
./ncrystal_ncmat2endf.py ncmat_file
```

with the following parameters:

```
positional arguments:
  input                 NCMAT file to convert

options:
  -h, --help            show this help message and exit
  -n NAME, --name NAME  Name of the compound in the NCMAT file (default: None)
  -t TEMPERATURES [TEMPERATURES ...], --temperatures TEMPERATURES [TEMPERATURES ...]
                        Temperatures to process the NCMAT file (default: [293.6])
  -e {greater,scaled,mixed}, --elastic_mode {greater,scaled,mixed}
                        Approximation used for the elastic component (default: scaled)
  -v VERBOSITY, --verbosity VERBOSITY
                        Controls how verbose should be the output (default: 1)
  -l {1,2,3,4,5}, --luxury {1,2,3,4,5}
                        Set the NCrystal vdoslux parameter used to generate the library. 3 is normal, 4 is fine and 5 is very fine. (default: 3)
  -g, --gif             Include the generalized information file (MF=7/MT=451) (default: False)
  -i, --isotopic        Expand each scatterer element into its isotopes (default: False)
  -m MATS, --mats MATS  JSON dictionary containing material number assignement for each element, e.g. '{"C":37, "H": 38}' (default: None)
  --alab ALAB           Set the ALAB parameter in MF1/MT451 (default: MyLAB)
  --auth AUTH           Set the AUTH parameter in MF1/MT451 (default: NCrystal)
  --libname LIBNAME     Set the LIBNAME parameter in MF1/MT451 (default: MyLib)
  --nlib NLIB           Set the NLIB parameter in MF1/MT451 (default: 0)
  --smin SMIN           Set the minimum value of S(alpha, beta) stored in MF7/MT4 (default: 1e-100)
  --lasym {0,1,2,3}     Write symmetric S(a,b) table (default: 0)

```
Ths script produces one TSL ENDF file per element in the material. The filenames are written as `tsl_SYMBOL_in_NAME.endf` where `NAME` is the
material name given with the `-n` parameter.

E.g., the following command creates a library for metallic aluminum at 350 K in the mixed elastic format in the file `tsl_Al.endf`, and assigns the ENDF material 53:

```
./ncrystal_ncmat2endf.py Al_sg225.ncmat -m '{"Al":53}' -n Al -t 350 -e mixed
```

The script can also be called from Python:


```
import ncrystal_ncmat2endf as ncmat2endf
endf_defaults = ncmat2endf.EndfParameters()
res = ncmat2endf.ncmat2endf('Al_sg225.ncmat', 'Al', endf_defaults, temperatures=[350], mat_numbers={'Al':53}, elastic_mode='mixed')
```

## Limitations
Currently the script is able to convert NCrystal polycrytalline materials defined with DOS or VDOS_Debye dynamic infos.

## Creating new evaluations

Many materials are currently available in the [NCrystal data library](https://github.com/mctools/ncrystal/wiki/Data-library). Additional materials can be created
using the [NCMATComposer](https://github.com/mctools/ncrystal-notebooks).

