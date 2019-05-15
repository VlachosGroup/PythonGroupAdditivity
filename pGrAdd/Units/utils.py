import numpy as np


def is_zero(val):
    # Compare scalar or array quantity to zero.
    if isinstance(val, np.ndarray):
        return (val == 0).all()
    else:
        return val == 0
