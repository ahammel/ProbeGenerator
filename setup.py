from distutils.core import setup
import sys

if sys.version < '3.2':
    print("probe_genertor requires Python v3.2 or later")
    sys.exit(1)

setup(name='ProbeGenerator',
      version='0.2.2',
      description='Automatic fusion-event probe maker',
      author='Alex Hammel',
      author_email='ahammel@bcgsc.ca',
      packages=['probe_generator'],
      requires=['docopt (>=0.6.1)',])
