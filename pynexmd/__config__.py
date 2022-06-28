import os, sys


DEBUG = False

MAX_MEMORY = int(os.environ.get('PYNEXMD_MAX_MEMORY', 4000)) # MB
TMPDIR = os.environ.get('TMPDIR', '.')
TMPDIR = os.environ.get('PYNEXMD_TMPDIR', TMPDIR)

VERBOSE = 3  # default logger level (logger.NOTE)
UNIT = 'angstrom'

#
# Loading pynexmd_conf.py and overwriting above parameters
#
for conf_file in (os.environ.get('PYNEXMD_CONFIG_FILE', None),
                  os.path.join(os.path.abspath('.'), '.pynexmd_conf.py'),
                  os.path.join(os.environ.get('HOME', '.'), '.pynexmd_conf.py')):
    if conf_file is not None and os.path.isfile(conf_file):
        break
else:
    conf_file = None

if conf_file is not None:
    if sys.version_info < (3,0):
        import imp
        imp.load_source('pynexmd.__config__', conf_file)
        del(imp)
    else:
        from importlib import machinery
        machinery.SourceFileLoader('pynexmd.__config__', conf_file).load_module()
        del(machinery)
del(os, sys)

