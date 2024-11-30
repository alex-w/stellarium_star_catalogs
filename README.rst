Codes to generate star catalogs for Stellarium
======================================================

**Work-in-progress**

This repository contains codes to generate experimental star catalogs for the popular open source planetarium software `Stellarium`_.
This experimental star catalogs are not meant to be used for astrophysical research, but rather aiming to provide a realistic 
sky for the Stellarium planetarium software.

The star catalogs are generated from the Hipparcos and Gaia DR3 catalogs.

.. _Stellarium: http://https://stellarium.org/

Requirements
------------

The code here are primarily written in Python using `Python>=3.12`. Run the following command to install the required Python packages:

```bash
pip install -r requirements.txt
```

You also need to install `embree`_ for high performance ray tracing.

.. _embree: https://www.embree.org/

Getting Data
------------

Most of the data are obtained from SIMBAD for Hipparcos stars and Gaia Archive for Gaia DR3 stars.

In order to generate synethic photometry (mainly `B-V` color for visual purposes) from Gaia DR3 continuous spectra, you need to download a large amount of data (~15 TB free disk space). The procudure is as follows:

#. | Download the data (DR3: 3362 GB) from this link - http://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/
   | To do this programmatically, you can run `wget -P ./ --no-clobber --no-verbose --no-parent --recursive --level=1 --no-directories http://cdn.gea.esac.esa.int/Gaia/gdr3/Spectroscopy/xp_continuous_mean_spectrum/`
#. | Unzip the data (DR3: 8204 GB)
   | To do this programmatically with 8 CPU cores, you can run `ls *.gz | xargs -n 1 -P 8 -I {} 'gunzip -c {} > ./$(basename {} .gz)'`
#. | Run GaiaXPy to synethize JKC photometry from the continuous spectra
   | To do this programmatically with 8 CPU cores, you can run `ls ./XpContinuousMeanSpectrum_2*.csv | xargs -n 1 -P 8 -I {} sh -c 'python gaiaxpy_phot.py {} ./jkc_photometry/tables'`

Scripts and Notebooks
----------------------

Here are the scripts and notebooks in this repository, you should follow the order to run the scripts and notebooks.

#. | `simbad_query_hip.py`_ - This script query SIMBAD for all possible Hipparcos IDs and save as a astropy Table.
#. | `Parse_HIP_binary.ipynb`_ - This notebook deal with issues like binary and mismatched HIP IDs with Gaia

.. _simbad_query_hip.py: simbad_query_hip.py
.. _Parse_HIP_binary.ipynb: Parse_HIP_binary.ipynb

Acknowledgement
----------------

If you are using this code to generate the star catalogs, please acknowledge the author and the following people/organizations:


Author
-------------
-  | **Henry Leung** - henrysky_
   | University of Toronto
   | Contact Henry: henryskyleung [at] gmail.com

License
-------

`GNU General Public License v3.0 <LICENSE>`_

The work done here is intended to be use in free and open source software (Stellarium and other open source softwares are welcome to adapt this code to generate star catalogs).