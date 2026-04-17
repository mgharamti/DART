
===============================
CMEMS SSH Along-Track Converter
===============================

.. contents:: 
   :depth: 3
   :local:

Overview
--------
This converter ingests **Copernicus Marine Environment Monitoring Service (CMEMS)**
along-track sea surface height (SSH) observations and converts them into a DART
``obs_seq`` file.

The supported product is:

- *Global Ocean Along Track L3 Sea Surface Heights NRT*
  (Product ID: ``SEALEVEL_GLO_PHY_L3_NRT_008_044``)

The converter reads CSV files exported from CMEMS and extracts **filtered sea level
anomaly (SLA)** observations for assimilation into ocean models such as ROMS.

The resulting observations are assigned the DART kind: ``SATELLITE_SSH``

Data Description
----------------
The converter uses the following variable: ``SLA_FILTERED``

This variable represents **low-pass filtered sea level anomaly**, including corrections for 
dynamic atmospheric correction (DAC), ocean tides, internal tides, and long-wavelength errors.
Compared to ``SLA_UNFILTERED``, this field is smoother, less noisy, and more 
representative of mesoscale variability. This makes it more appropriate for 
comparison with model sea surface height (e.g., ROMS ``zeta``).

Observation Error
-----------------
The observation error standard deviation is parameterized as a function of
**local bathymetry**:

::

   sigma = sigma_min + (sigma_max - sigma_min) * exp(-h / h0)

where:

- ``sigma_min``: minimum error (deep ocean)
- ``sigma_max``: maximum error (shallow water)
- ``h``: bathymetry depth (meters)
- ``h0``: transition depth scale

This formulation reflects increased uncertainty in shallow/coastal regions due to 
land contamination, unresolved tides (barotropic and internal), wetting/drying 
effects, stronger small-scale variability, and interpolation mismatch.

If bathymetry cannot be determined or is invalid:

::

   sigma = sigma_max

Bathymetry
^^^^^^^^^^
Bathymetry is read from a model grid file (e.g., ROMS restart or grid file).
The converter supports two modes: 

**1. ROMS fast lookup (default)**

- Uses 1D slices of longitude and latitude:
- Finds nearest indices and searches a small local window (5x5 cells)

This approach is efficient and suitable for structured ROMS grids.

**2. Fallback brute-force search**

- Used if the grid is not monotonic
- Searches the entire grid

Both methods ignore land points using the grid mask.

Time Handling
-------------
Observation times are read from ISO-formatted strings, for example:

::

   2026-04-16T08:43:54.000Z

These are converted to DART ``time_type`` using the Gregorian calendar.

Usage
-----
1. Prepare a list of input CSV files and place them in a text file, for instance ``ssh_list.txt`` 

2. Configure the namelist in ``input.nml``:

::

   &cmems_ssh_to_obs_nml
      file_list         = 'ssh_file_list.txt'
      ocean_in          = 'roms_restart.nc'
      file_out          = 'obs_seq.ssh'
      avg_obs_per_file  = 500000
      obs_error_min     = 0.04
      obs_error_max     = 0.08
      transition_depth  = 500.0
      debug             = .true.
   /

+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| Namelist item        | Type                  | Default Values            | Description                                                  |
+======================+=======================+===========================+==============================================================+
| ``file_list``        | character(len=256)    | ``''`` (empty string)     | Path to a text file containing a list of input CSV files,    |
|                      |                       |                           | one filename per line.                                       |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``ocean_in``         | character(len=256)    | ``'roms_restart.nc'``     | Path to ocean model file used to read bathymetry, grid       |
|                      |                       |                           | coordinates, and land mask.                                  |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``file_out``         | character(len=256)    | ``'obs_seq.ssh'``         | Name of the DART observation sequence output file. If the    |
|                      |                       |                           | file already exists, it will be replaced.                    |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``avg_obs_per_file`` | integer               | ``500000``                | Estimated average number of observations per input file.     |
|                      |                       |                           | Used to pre-allocate sequence memory efficiently.            |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``obs_error_min``    | real(r8)              | ``0.04``                  | Minimum observation error standard deviation (meters),       |
|                      |                       |                           | typically applied in deep ocean regions.                     |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``obs_error_max``    | real(r8)              | ``0.08``                  | Maximum observation error standard deviation (meters),       |
|                      |                       |                           | typically applied in shallow/coastal regions.                |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``transition_depth`` | real(r8)              | ``500.0``                 | Bathymetric transition depth (meters) controlling the        |
|                      |                       |                           | exponential error scaling.                                   |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+
| ``debug``            | logical               | ``.true.``                | If true, prints diagnostic output including sample           |
|                      |                       |                           | observations and processing information.                     |
+----------------------+-----------------------+---------------------------+--------------------------------------------------------------+

3. Run the converter: ``./cmems_ssh_to_obs``

4. Output: The resulting observation sequence file is stored in ``obs_seq.ssh``

When ``debug`` is turned on, the progress of the converter can be monitored 
which would look similar to the following:

::

  Input file: #2 ../data/obs_ssh_2.csv

   CSV Header Fields (10 columns):
     1  parameter
     2  platformId
     3  platformType
     4  time
     5  longitude
     6  latitude
     7  depth
     8  pressure
     9  value
    10  valueQc

    * obs #  1, lat: -14.2924, lon: 127.8070, SSH: 0.2940, QC: 0, date: 2026 Apr 14 19:26:36
    * obs #  2, lat: -14.2351, lon: 127.8010, SSH: 0.2820, QC: 0, date: 2026 Apr 14 19:26:37
    * obs #  3, lat: -14.1778, lon: 127.7949, SSH: 0.2690, QC: 0, date: 2026 Apr 14 19:26:38
    * obs #  4, lat: -14.1206, lon: 127.7888, SSH: 0.2560, QC: 0, date: 2026 Apr 14 19:26:39
    * obs #  5, lat: -14.0633, lon: 127.7827, SSH: 0.2420, QC: 0, date: 2026 Apr 14 19:26:40
    * obs #  6, lat: -14.0060, lon: 127.7766, SSH: 0.2280, QC: 0, date: 2026 Apr 14 19:26:41
    * obs #  7, lat: -13.9487, lon: 127.7705, SSH: 0.2130, QC: 0, date: 2026 Apr 14 19:26:42
    * obs #  8, lat: -13.8915, lon: 127.7645, SSH: 0.2000, QC: 0, date: 2026 Apr 14 19:26:43
    * obs #  9, lat: -13.8342, lon: 127.7584, SSH: 0.1910, QC: 0, date: 2026 Apr 14 19:26:43
    * obs # 10, lat: -13.7769, lon: 127.7523, SSH: 0.1850, QC: 0, date: 2026 Apr 14 19:26:44

 > Ready to write 3042 observations:
   write_obs_seq  opening formatted observation sequence file "obs_seq.ssh"
   write_obs_seq  closed observation sequence file "obs_seq.ssh"
   cmems_ssh_to_obs Finished successfully.

References
----------
- Copernicus Marine Service:
  https://marine.copernicus.eu

- Product documentation:
  https://data.marine.copernicus.eu/product/SEALEVEL_GLO_PHY_L3_NRT_008_044
