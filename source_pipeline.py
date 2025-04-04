"""Compute subproducts from data cubes.

Not all the steps run were published.
"""
from configparser import ConfigParser
from pathlib import Path
import sys

from astro_source.source import Source
from line_little_helper.moving_moments import moving_moments
from line_little_helper.line_peak_map import line_peak_map
from line_little_helper.molecule import NoTransitionError
from line_little_helper.pvmap_extractor import pvmap_extractor
from line_little_helper.subcube_extractor import subcube_extractor
from line_little_helper.symmetric_moments import symmetric_moments
from line_little_helper.utils import normalize_qns
from tile_plotter.plotter import plotter


from common_paths import MOL_DIR, CONFIGS, RESULTS, FIGURES

# Recommended line lists
LINE_LISTS = {
#    '(13)CO': 'JPL',
#    'H2CO': 'JPL',
#    'SiO': 'JPL',
#    'Hα': 'Recomb',
#    'CH3CN': 'JPL',
#    '(13)CH3CN': 'CDMS',
    'CH3OH': 'CDMS',
#    '(13)CH3OH': 'CDMS',
#    'CH3CHO': 'SLAIM',
#    'CH3OCHO': 'SLAIM',
#    'HNCO': 'CDMS',
    'SO': 'JPL',
#    'SO2': 'CDMS',
#    'HC3N': 'CDMS',
#    'NH2CHO': 'CDMS',
}

# Molecules and transitions of interest
LINE_TRANSITIONS = {
    'CH3OH': [#'4(2,3)-5(1,4)A,vt=0',       #spw0
              #'5(4,2)-6(3,3)E,vt=0',
              '18(3,15)-17(4,14)A,vt=0',
              #'10(2,9)-9(3,6)A,vt=0',      #spw1
              #'10(2,8)-9(3,7)A,vt=0',
              #'18(3,16)-17(4,13)A,vt=0',
              #'4(-2,3)-3(-1,2)E,vt=0',      #spw2
              #'5(-1,4)-4(-2,3)E,vt=0',
              #'20(-1,19)-20(-0,20)E,vt=0',
              #'8(-0,8)-7(-1,6)E,vt=0',     #spw3
              #'23(-5,18)-22(-6,17)E,vt=0',
              #'25(-3,23)-24(-4,20)E,vt=0',
              #'6(1,5)-7(2,6)--,vt=1',      # other
              #'13(3,10)-14(4,10)A,vt=2',
              #'34(-13,22)-33(-11,22)E,vt=0-vt=1',
              ],
    #'(13)CH3OH': [#'10(2,8)-9(3,7)++',
    #              '5(1,5)-4(1,4)++'],
    #'CH3CN': [#'12(0)-11(0)',
    #          #'12(1)-11(1)',
    #          '12(2)-11(2)',
    #          '12(3)-11(3)',
    #          '12(4)-11(4)',
    #          '12(8)-11(8)'],
    #'(13)CH3CN': ['13(3)-12(3)',
    #              '13(4)-12(4)'],
    #'CH3CHO': ['2(2,0)-3(1,3)A++,vt=2'],
    #'HNCO': ['28(1,28)-29(0,29)',
    #         '10(3,8)-9(3,7)',
    #         '10(0,10)-9(0,9)'],
    #'CH3OCHO': ['20(0,20)-19(0,19)A',
    #            '32(9,24)-32(8,25)E',
    #            '48(14,35)-47(15,32)E'],
    #'SiO': ['5-4'],
    'SO': ['6(5)-5(4)'],
    #'SO2': ['59(14,46)-60(13,47)',
    #        '28(3,25)-28(2,26)',
    #        '99(9,91)-98(10,88)',
    #        '94(21,73)-95(20,76)'],
    #'HC3N': ['J=24-23,l=0',
    #         'J=24-23,l=2e',
    #         'J=24-23,l=1f'],
    #'NH2CHO': ['11(2,10)-10(2,9)'],
}

# Saved molecules
SAVED_MOLS = {
    'CH3OH': MOL_DIR / 'ch3oh.json',
    #'CH3CN': MOL_DIR / 'ch3cn.json',
    #'(13)CH3OH': MOL_DIR / '13ch3oh.json',
    #'(13)CH3CN': MOL_DIR / '13ch3cn.json',
    #'SiO': MOL_DIR / 'sio.json',
    'SO': MOL_DIR / 'so.json',
}

def search_molecule(source, molecule, array):
    data = []
    for key in src.get_data_sections():
        if 'molecules' not in src.config[key] or array not in key:
            continue
        molecules = src.config[key]['molecules'].split()
        if molecule in molecules:
            data.append(src.config[key])

    return data

def gen_plot_config(source, moments, template, outfile, mol=None, qns=None):
    # Open template
    if outfile.exists():
        return
    config = ConfigParser()
    config.read(template)

    # Moment sections
    for i in range(3):
        section = f'moment{i}'
        position = source.position.to_string(style='hmsdms')
        position = f'{position} {source.position.frame.name}'
        config[section]['moment'] = str(moments[i])
        config[section]['center'] = position
        config[section]['label'] = f'Moment {i}'
        if i == 0 and mol is not None:
            title = f'{mol}' if qns is None else f'{mol} ({qns})'
            config[section]['title'] = title
        if i == 1:
            # Shift data is additive
            vlsr = f'{-source.vlsr.value}'
            vlsr += ' ' + f'{source.vlsr.unit}'.replace(' ', '')
            config[section]['shift_data'] = vlsr

    # Continuum
    config['continuum']['contour'] = source.config['b6_c8c9_continuum']['file']

    # Write
    with open(outfile, 'w') as configfile:
        config.write(configfile)

def crop_line(source,
              outdir,
              configs,
              figures,
              array,
              molecules=['CH3OH'],
              qns_mol = LINE_TRANSITIONS,
              half_width=30):
    for mol in molecules:
        configs_with_mol = search_molecule(source, mol, array)
        processed = []
        norm_mol = mol.replace('(', '').replace(')', '')
        for qns in qns_mol[mol]:
            if qns in processed:
                continue
            norm_qns = normalize_qns(qns)

            for src_cfg in configs_with_mol:
                cube = src_cfg['file']
                
                # Flags
                flags = ['--vlsr', f'{source.vlsr.value}',
                         f'{source.vlsr.unit}'.replace(' ', ''),
                         '--line_lists', line_lists.get(mol, 'CDMS'),
                         '--molecule', mol,
                         '--qns', qns,
                         '--common_beam',
                         '--put_linefreq',
                         '--win_halfwidth', f'{half_width}',
                         '--spectral_axis', 'velocity',
                        ]
                if 'rms' in src_cfg:
                    flags += ['--rms'] + src_cfg['rms'].split()

                # Is there a mask?
                moldir = outdir / norm_mol
                mask = moldir / 'source_mask.fits'
                if mask.is_file():
                    flags += ['--mask', f'{mask}', '--shrink']
                    masked = '_masked_shrink'
                else:
                    masked = ''
                
                # Compute moments
                try:
                    print(cube)
                    name = f'{norm_mol}_{norm_qns}_nchans{half_width*2}{masked}.fits'
                    out_cube = moldir / name
                    moldir.mkdir(parents=True, exist_ok=True)
                    subcube_extractor(flags + [cube, f'{out_cube}'])
                    processed.append(qns)
                except NoTransitionError:
                    print(f'{mol} ({qns}): not in cube {cube}')
                    continue

def moments(source,
            outdir,
            configs,
            figures,
            array,
            molecules=['CH3OH'],
            qns_mol=LINE_TRANSITIONS,
            half_width=20,
            nsigma=3):
    plot_template = configs / 'templates' / 'moment_maps.cfg'
    for mol in molecules:
        # Find what molecules are in the source
        configs_with_mol = search_molecule(source, mol, array)
        processed = []
        norm_mol = mol.replace('(', '').replace(')', '')
        for qns in qns_mol[mol]:
            if qns in processed:
                continue

            for src_cfg in configs_with_mol:
                cube = src_cfg['file']
                # Plot config
                suffix = f'{src_cfg.name}_width{half_width*2}_nsig{nsigma}'
                norm_qns = normalize_qns(qns)
                #cfg = f'{source.name}_{array}_{norm_mol}_{norm_qns}.cfg'
                cfg = f'{source.name}_{norm_mol}_{norm_qns}_{suffix}.cfg'
                cfg = configs / 'plots/processing' / cfg
                cfg.parent.mkdir(parents=True, exist_ok=True)
                if cfg.exists():
                    continue
                
                # Moment flags
                flags = ['--vlsr', f'{source.vlsr.value}',
                         f'{source.vlsr.unit}'.replace(' ', ''),
                         '--line_lists', LINE_LISTS.get(mol, 'CDMS'),
                         '--molecule', mol,
                         '--qns', qns,
                         '--nsigma', f'{nsigma}',
                         '--win_halfwidth', f'{half_width}',
                        ]
                if mol in SAVED_MOLS:
                    flags += ['--restore_molecule', f'{SAVED_MOLS[mol]}']
                if 'rms' in src_cfg:
                    flags += ['--rms'] + src_cfg['rms'].split()
                
                # Compute moments
                try:
                    name = f'{norm_mol}_{norm_qns}_{suffix}'
                    moldir = outdir / norm_mol
                    moldir.mkdir(parents=True, exist_ok=True)
                    filenames = symmetric_moments(
                        flags + [cube, str(moldir / name), '0', '1', '2'])
                    processed.append(qns)
                except NoTransitionError:
                    print(f'{mol} ({qns}): not in cube {cube}')
                    continue

                # Generate plot config
                gen_plot_config(source, filenames, plot_template, cfg, mol=mol,
                                qns=qns)

                # Plot
                plotname = figures / source.name / array
                plotname.mkdir(parents=True, exist_ok=True)
                plotname = plotname / f'{mol}_{norm_qns}_{suffix}.png'
                plotter([f'{cfg}', f'{plotname}'])

def pv_maps(source, outdir, configs, figures, array):
    # Streamers pv
    streams = configs / f'pvmaps/{source.name}_streams.cfg'
    if streams.is_file():
        pv_streams(streams, source, outdir, figures, array)

    # Rotation pvs
    configpv_rot = configs / 'pvmaps' 
    for rotation in configpv_rot.glob(f'{source.name}_rotation*.cfg'):
        try:
            print(rotation)
            pv_rotation(rotation, source, outdir, figures, array)
        except OSError:
            print('Error')
            continue

def pv_extractor(pvconfig, source, outdir, recenter=False):
    results = outdir / f'{pvconfig.stem}.fits'
    flags = ['--pvconfig', f'{pvconfig}', '--output', f'{results}',
             '--estimate_error', '--source', f'{source.config_file}',
             '--common_beam']
    if recenter:
        flags += ['--recenter']
    return pvmap_extractor(flags)
    #return results.parent.glob(f'{pvconfig.stem}*.fits')
    
def pv_streams(pvconfig, source, outdir, figures, array):
    # Extract pv maps
    pvout = outdir / 'pvmaps'
    pvmaps = list(pvout.glob('*stream*.fits'))
    if len(pvmaps) == 0:
        pvmaps = pv_extractor(pvconfig, source, pvout)

def pv_rotation(pvconfig, source, outdir, figures, array):
    pvout = outdir / 'pvmaps'
    pvmaps = pv_extractor(pvconfig, source, pvout, recenter=True)

def split_moments(source, outdir, configs, figures, array,
                  molecules=['SO'],
                  chansep=9,
                  chanwidth=5):
    # Calculate split moments
    for mol in molecules:
        configs_with_mol = search_molecule(source, mol, array)
        for config in configs_with_mol:
            cube = config['file']
            flags = ['--vlsr', f'{source.vlsr.value}',
                     f'{source.vlsr.unit}'.replace(' ', ''),
                     '--line_lists', LINE_LISTS.get(mol, 'CDMS'),
                     '--savemasks',
                     '--split', f'{chansep}', f'{chanwidth}',
                     '--molecule', mol]
            if mol in SAVED_MOLS:
                flags += ['--restore_molecule', f'{SAVED_MOLS[mol]}']
            if 'rms' in config:
                flags += ['--rms'] + config['rms'].split()
            dir_suff = f'{mol}_split{chansep}_{chanwidth}_{config.name}'
            moving_moments(flags + ['80', f'{outdir / dir_suff}', cube])

def peak_maps(source,
              outdir,
              configs,
              figures,
              array,
              molecules=['CH3OH'],
              qns_mol={'CH3OH':['18(3,15)-17(4,14)A,vt=0']},
              half_width=10):
    for mol in molecules:
        configs_with_mol = search_molecule(source, mol, array)
        processed = []
        norm_mol = mol.replace('(', '').replace(')', '')
        for qns in qns_mol[mol]:
            if qns in processed:
                continue

            for src_cfg in configs_with_mol:
                cube = src_cfg['file']
                # Plot config
                norm_qns = normalize_qns(qns)
                
                # Moment flags
                flags = ['--vlsr', f'{source.vlsr.value}',
                         f'{source.vlsr.unit}'.replace(' ', ''),
                         '--line_lists', line_lists.get(mol, 'CDMS'),
                         '--molecule', mol,
                         '--qns', qns,
                         '--nsigma', '5',
                         '--win_halfwidth', f'{half_width}',
                         '--moments', '0', '1',
                         '--nlinewidth', '2',
                         '--use_dask',
                         '--common_beam',
                        ]
                if 'rms' in src_cfg:
                    flags += ['--rms'] + src_cfg['rms'].split()
                
                # Compute moments
                try:
                    name = f'{norm_mol}_{norm_qns}'
                    moldir = outdir / norm_mol
                    moldir.mkdir(parents=True, exist_ok=True)
                    filenames = line_peak_map(flags + [cube, f'{moldir}'])
                    processed.append(qns)
                except NoTransitionError:
                    print(f'{mol} ({qns}): not in cube {cube}')
                    continue

if __name__ == '__main__':
    ## Constants
    figures = Path('/data/share/dihca2/combined_projects/figures')

    # Steps
    steps = {
        1: moments,
        2: split_moments,
        3: pv_maps,
        4: crop_line,
        5: peak_maps,
    }
    skip = [1, 2, 4, 5]
    array = 'c8c9'

    # Read sources from command line
    sources = ['G336.01-0.82']

    # Iterate over source config files
    iterover = (CONFIGS / f'sources/{source}.cfg' for source in sources)
    for config in iterover:
        # Open source
        print(config)
        src = Source(config_file=config)
        outdir = RESULTS / src.name / array

        # Run steps
        for n, func in steps.items():
            if n in skip:
                continue
            print(f'Step {n}')
            func(src, outdir, CONFIGS, FIGURES, array)
