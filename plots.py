from pathlib import Path

from tile_plotter.plotter import plotter

#from common_paths import *

plot_type = 'papers'
#plot_type = 'presentations'

if __name__ == '__main__':
    # Constants
    results = Path('/data/share/dihca2/combined_projects/results')
    configs = Path('/data/share/dihca2/combined_projects/scripts/configs')
    figures = Path('/data/share/dihca2/combined_projects/figures')
    config_dir = configs / 'plots' / plot_type

    # Steps
    steps = {
        1: 'continuum',
        2: 'moments',
        3: 'split_moments',
        4: 'pv_maps',
        5: 'streamers',
        6: 'composite',
    }
    skip = [3, 2, 4, 5, 6]

    # Read sources from command line
    sources = ['G336.01-0.82']

    # Flags
    flags = []
    if plot_type == 'papers':
        flags = ['--pdf']

    # Iterate over sources
    for source in sources:
        # Iterate over steps
        for key, val in steps.items():
            if key in skip:
                continue

            config = config_dir / f'{source}_{val}.cfg'
            if config.exists():
                print(f'Plotting {config}')
                plotname = figures / source / plot_type / f'{val}.png'
                plotter([f'{config}', f'{plotname}'] + flags)
            elif itercfgs := config_dir.glob(f'{source}_{val}_*.cfg'):
                for config in itercfgs:
                    print(f'Plotting {config}')
                    name = config.stem.split('_')
                    ind = name.index(val.split('_')[0])
                    name = '_'.join(name[ind:])
                    plotname = figures / source / plot_type / f'{name}.png'
                    plotter([f'{config}', f'{plotname}'] + flags)
            else:
                print(f'Skipping: {config}')
                continue
