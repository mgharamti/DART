.. _adding-your-observations-to-dart:

Adding your observations to DART
================================

First, you should understand that DART already supports a tremendous variety of
observations. To fully support an observation means that the observation can be
converted from its native format to the DART observation sequence format and
that the observation forward operator is already implemented. Keep in mind that
forward operators are not specific to any one model.

The next sections describe how observations are handled in DART.
Following those are sections grouped under "Observation Converters"
describing types of observations; their sources and methods for 
formatting them for use by DART.

DART can use both real and "synthetic" observations.
Synthetic observations can be created manually using 
:ref:`create_obs_sequence`
or extracted from a numerical model using perfect_model_obs 
as described in *Identity observations*, below.


Real observations
-----------------

The observation converters are in the *observations/obs_converter* directory,
and are documented in :ref:`available_observation_converters`. 

The forward operators are functionally or logically grouped into Fortran modules
in the *observations/forward_operator* directory. DART employs a ‘contractual’
style of programming in that the forward operator requests information from the
model, and if the model cannot provide it, the forward operator may request
different information in an attempt to collect the information needed to apply
the operator. If the model cannot provide any of the required information, the
forward operator fails, the DART QC for that observation is set to the
appropriate value, and the program continues.


Identity observations from perfect_model_obs
--------------------------------------------

An identity observation is a special type of observation, 
which samples a model state directly.
It is a copy of a chosen state variable. 
That is, the forward operator H(x) is the identity matrix (hence the name). 
The value of the observation is taken from the model state
at a chosen grid point, which is the location of the observation.
The observation error variance is specified for later use
and an observation error is added to the value; 
typically a random draw from an distribution sized according to the error variance.
Identity observations do not get listed as a type in the header of an 
observation sequence file, they are denoted in a given observation by a type
of -x where x is the index in the DART state vector that the observation 
corresponds to.
This means that identity observations from one model can't (usually) be used
in an assimilation with another model.
They are useful for OSSEs and model interface development.
The process for creating identity observations is outlined in
:ref:`creating-obs-seq-synthetic` and 
:ref:`synthetic_observations`.
