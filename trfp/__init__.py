"""Module for doing trolley-fixed probe correlations."""

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
from trfp.geometry import STATION_RING_PHI
from trfp.geometry import STATION_PROBE_ID
from trfp.geometry import PLUNGING_PROBE_CALIBRATIONS

import trfp.matrices
from trfp.matrices import THETA_FP_4
from trfp.matrices import THETA_FP_4_ST37_ST39
from trfp.matrices import THETA_FP_4_ST41
from trfp.matrices import THETA_FP_6
from trfp.matrices import THETA_TR
from trfp.matrices import J_6_PROBE
from trfp.matrices import J_6_PROBE_OFFSET
from trfp.matrices import J_4_PROBE
from trfp.matrices import J_4_PROBE_ST37_ST39
from trfp.matrices import J_4_PROBE_ST41

import trfp.runs
from trfp.runs import Run

import trfp.analysis
from trfp.analysis import remove_trolley_effect
from trfp.analysis import trolley_run_station_average

import trfp.uncertainty