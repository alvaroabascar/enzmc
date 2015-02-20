from distutils.core import setup, Extension

monty = Extension('monty', ['monty.c'])

setup(name = 'monty',
      version = '0.0',
      description = 'monty module',
      ext_modules = [monty])
