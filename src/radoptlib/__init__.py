"""RadOptLIB: a library of radiotherapy treatment-planning optimization problem instances."""

from .convert import to_mps
from .io import Problem, data_dir, download, info, list_instances, load, load_solution, select
from .kkt_metrics import verify_solution

__version__ = "1.0.0"

__all__ = ['Problem', 'data_dir', 'download', 'info', 'list_instances', 'load', 'load_solution', 'select',
           'to_mps', 'verify_solution']