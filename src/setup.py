from distutils.core import setup, Extension

sources = ['monty.c', 'montecarlo/montecarlo.c', 'nlr/lvmrq.c',
           'random/random.c', 'misc/mathlib.c', 'misc/matrix.c',
           'lineq/gaussjbs.c', 'models/models.c']
monty = Extension('monty',
                  include_dirs = ['/usr/local/include', './include/'],
                  sources = sources)

setup(name = 'monty',
      version = '0.0',
      description = 'monty module',
      ext_modules = [monty])
