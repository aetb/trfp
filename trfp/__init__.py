"""Module for doing trolley-fixed probe correlations."""

## Geometry sub-module

import trfp.geometry
from trfp.geometry import FP4_X
from trfp.geometry import FP4_X_ST37_ST39
from trfp.geometry import FP4_X_ST41
from trfp.geometry import FP4_Y
from trfp.geometry import FP6_X
from trfp.geometry import FP6_X_OFFSET
from trfp.geometry import FP6_Y
from trfp.geometry import TR_X
from trfp.geometry import TR_Y
from trfp.geometry import STATION_BARCODE_PHI
from trfp.geometry import STATION_BARCODE_EDGES
from trfp.geometry import STATION_BARCODE_PHI_6
from trfp.geometry import STATION_BARCODE_EDGES_6
from trfp.geometry import STATION_RING_PHI
from trfp.geometry import STATION_PROBE_ID
from trfp.geometry import STATION_PROBE_NUMBER
from trfp.geometry import PLUNGING_PROBE_CALIBRATIONS

## Matrices sub-module

import trfp.matrices
from trfp.matrices import THETA_FP_4
from trfp.matrices import THETA_FP_4_ST37_ST39
from trfp.matrices import THETA_FP_4_ST41
from trfp.matrices import THETA_FP_6
from trfp.matrices import THETA_TR
from trfp.matrices import J_6_PROBE
from trfp.matrices import J_6_PROBE_OFFSET
from trfp.matrices import J_ST_5
from trfp.matrices import J_4_PROBE
from trfp.matrices import J_4_PROBE_ST37_ST39
from trfp.matrices import J_4_PROBE_ST41


## Runs sub-module
## This sub-module needs major renovations for v3.

import trfp.runs
from trfp.runs import Run

## Analysis sub-module
## This module should contain a lot of the helper_functions.

import trfp.analysis
from trfp.analysis import *

## Plotting sub-module
## Has a few useful plotting methods.

import trfp.plotting

## Configuration files
## Dictionaries that define the different runs of Run 1, 2, and 3.

import trfp.field_map_config_run1
import trfp.field_map_config_run2
import trfp.field_map_config_run3