# imports
# ----------------------------------------------------------------------------------------------------------------------
import sys
from pathlib import Path
# ----------------------------------------------------------------------------------------------------------------------


def log_out(path_log: (str, Path), log_name: str = f"{Path(__file__).stem}__log"):
    """
    Simple logging function for saving the stdout and stderr to a .log file.

    Parameters
    ----------
    path_log: str | Path
        Path for saving the .log file.
    log_name: str
        Name of the .log file.
    """
    path_log = Path(path_log)
    log_file = path_log / f"{log_name}.log"
    if not log_file.exists():
        log_file.parent.resolve().mkdir(exist_ok=True, parents=True)

    # create log-file and set the standard system out and error to log
    log = open(log_file, "a")
    sys.stdout = log
    sys.stderr = log

