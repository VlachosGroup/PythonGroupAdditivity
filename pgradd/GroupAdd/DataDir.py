import os
import sys

# What environment variable is used to override use of the bundled
# data directory?
_data_dir_envvar = 'pgradd_DATA_DIR'

# We'll cache the located data directory root in this variable:
_data_dir_cached = False


def get_data_dir():
    """Returns the directory containing the pgradd data library.

    Return the cached directory path if it has been set.  Otherwise:

    1. If the _data_dir_envvar is set in the runtime environment,
       cache and use it as the data dir

    2. This Python file exists somewhere under the 'pgradd' module
       directory; step back until we get a path that ends with '/pgradd'
       and tack '/data' onto the end of it.  Voila, our bundled data
       directory.

    The caveat to 2. is when this script is run from the GroupAdd
    directory for testing purposes; the __file__ will be just
    "DataDir.py" so it needs to be made absolute.
    """
    global _data_dir_cached, _data_dir_cached

    if not _data_dir_cached:
        base_path = os.getenv(_data_dir_envvar, False)
        if not base_path:
            base_path = os.path.abspath(os.path.dirname(__file__))
            while base_path != os.sep and not base_path.endswith(os.sep +
                                                                 'pgradd'):
                base_path = os.path.dirname(base_path)
            if base_path == os.sep:
                raise RuntimeError('DataDir.get_data_dir: unable to locate',
                                   'bundled data directory')
            base_path = os.path.join(base_path, 'data')
        if not os.path.isdir(base_path):
            raise RuntimeError('DataDir.get_data_dir: directory does not',
                               'exist: ' + base_path)
        _data_dir_cached = base_path
    return _data_dir_cached

#
# This code can be tested very simply by
#
#   python DataDir.py
#


if __name__ == '__main__':
    sys.stdout.write('pgradd data directory = {0:s}'.format(get_data_dir()))
