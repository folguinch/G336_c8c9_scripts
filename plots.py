"""Plot data.

The global parameter `plot_type` specify the subdirectory where the
configuration files for the plot type are available. The filename of the
configuration files follow the pattern `<source name>_<step>.cfg` or
`<source name>_<step>_<suffix>.cfg`, where the available `step` values are
listed in the parameter `steps`. Steps can be selected with with the `skip`
parameters.
"""
from pathlib import Path

from tile_plotter.plotter import plotter

from common_paths import CONFIGS, FIGURES

plot_type = 'papers'
#plot_type = 'presentations'

if __name__ == '__main__':
    # Constants
    config_dir = CONFIGS / 'plots' / plot_type

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
                plotname = FIGURES / source / plot_type / f'{val}.png'
                plotter([f'{config}', f'{plotname}'] + flags)
            elif itercfgs := config_dir.glob(f'{source}_{val}_*.cfg'):
                for config in itercfgs:
                    print(f'Plotting {config}')
                    name = config.stem.split('_')
                    ind = name.index(val.split('_')[0])
                    name = '_'.join(name[ind:])
                    plotname = FIGURES / source / plot_type / f'{name}.png'
                    plotter([f'{config}', f'{plotname}'] + flags)
            else:
                print(f'Skipping: {config}')
                continue
