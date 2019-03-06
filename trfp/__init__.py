"""Module for doing trolley-fixed probe correlations."""

from trfp.geometry import FP_4_X
from trfp.geometry import FP_4_Y
from trfp.geometry import FP_6_X
from trfp.geometry import FP_6_Y
from trfp.geometry import TR_X
from trfp.geometry import TR_Y
from trfp.geometry import STATION_BARCODE_PHI
from trfp.geometry import STATION_RING_PHI
from trfp.geometry import STATION_PROBE_ID

from trfp.matrices import THETA_FP_4
from trfp.matrices import THETA_FP_6
from trfp.matrices import THETA_TR

from trfp.runs import TrolleyRun
from trfp.runs import FixedProbeRun

from trfp.analysis import remove_trolley_effect