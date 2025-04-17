#!/bin/python3
"""Update to FERIA model based on new centrifugal radius.

Centrifugal barrier = 0.5 * centrifugal radius

Before runing make in FERIA, modify the lbNpix and lbNvel to 8.
Change the fname_exec parameter accordingly to your implementation of FERIA.
"""
from itertools import product
from pathlib import Path
from string import Template
import sys
import subprocess
import math

import numpy as np

from common_paths import RESULTS

# Change path to FERIA according to your system
fname_exec = '~/clones/feria/feria'
basedir = RESULTS / 'G336.01-0.82/c8c9/CH3OH'
fname_cubein = Path('./feria.in')

# This is to compare with previous results, so we do not need to update the
# position
source = {'object': 'G336',
          'ra': '16h35m09.261s',
          'dec': '-48d46m47.66s',
          'vsys': '-47.2'}

obs_pars = {'line': 'CH3OH',
            'restfreq': '233.795666',
            'pixsize': '0.004',
            'velres': '0.1'}

phys_pars = {'distance': ['3100'], # pc
             'mass': ['10'], # msun
             'rcb': ['250'], # au
             'incl': ['65.'],
             #'pa': ['125'],
             'pa': ['-55'],
             #'rot': ['1'],
             'rot': ['-1'],
             'rout': ['800.'],
             'rin': ['CB'],
             'ireheight': ['0.'],
             'ireflare': ['30'],
             'irenprof': ['-1.5'],
             'iretprof': ['-0.4'],
             'kepheight': ['0.5'],
             'kepflare': ['30'],
             'kepnprof': ['-1.5'],
             'keptprof': ['-0.4'],
             'cbdens': ['1e-2'],
             'cbtemp': ['50'],
             'lw': ['1.0'],
             'bmaj': ['0.0692'],
             'bmin': ['0.0517'],
             'bpa': ['-29.07']
             }

pv_pars = {'pvpa': ['125'],
           'pvra': ['0.0'],
           'pvdec': ['0.0']}



####################

def make_cube(pars, output, template=Path('./feria_template.in')):
    if pars['rin'] == 'CB':
        pars['rin'] = pars['rcb']
    template = Template(template.read_text())
    output.write_text(template.substitute(**pars))
    subprocess.call(f'{fname_exec} < {output}', shell = True)

def normalize_name(name, **pars):
    return f'{name}_' + '_'.join(f'{k}{v}' for k, v in pars.items()) + '.fits'


ncube = np.prod([len(val) for val in phys_pars.values()])
npv = np.prod([len(val) for val in pv_pars.values()])

for i, physvals in enumerate(product(*phys_pars.values())):
    print(f'<<<<< Make a Cube ({i+1} / {ncube}) >>>>>')
    pars = dict(zip(phys_pars.keys(), physvals))
    cubename = normalize_name(source['object'], **pars)
    pars.update(source)
    pars.update(obs_pars)

    # Separate by disk PA
    modeldir = basedir / f"feria_models_PA{pars['pa']}"
    modeldir.mkdir(exist_ok=True)
    cubename = modeldir / cubename
    if cubename.exists():
        print('Model already calculated, skipping')
        continue
    

    for j, pvvals in enumerate(product(*pv_pars.values())):
        print(f'<<< Make a PV ({j+1} / {npv}) >>>')
        pvdict = dict(zip(pv_pars.keys(), pvvals))
        #pvname = modeldir / normalize_name(cubename.stem, **pvdict)
        pars.update(pvdict)
        pars['cubename'] = f'{modeldir}/feria_model_updated.fits'

        make_cube(pars, fname_cubein)

    #    # Rename files
    #    print(f'Moving files to {modeldir}')
    #    pvfiles = list(basedir.glob('*PV*.fits'))
    #    if len(pvfiles) > 1:
    #        raise ValueError('Too many pv maps')
    #    pvfiles[0].rename(pvname)
    #    cubefiles = list(basedir.glob('*.fits'))
    #    if len(cubefiles) > 1:
    #        raise ValueError('Too many cubes')
    #    cubefiles[0].rename(cubename)

print('\n\007\007')
print(f'\n\n\n----------\n{sys.argv[0]} has been done.\n\n\n')


