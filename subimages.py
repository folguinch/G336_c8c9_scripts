"""Sub-images for online storage."""
from pathlib import Path
from casatasks import imsubimage, exportfits, imhead

outdir = Path('../data/cubes/subimages')
indir = Path('../data/cubes/')
images = [
    indir / 'G336.018-0.827_spw0_CH3OH_robust2_multiscale.image.fits',
    indir / 'G336.018-0.827_spw2_SO_robust2_multiscale.image.fits',
    indir / 'G336.018-0.827_spw3_SO_robust2_multiscale.image.fits']
box = '3100,3100,7700,7700'
outdir.mkdir(parents=True, exist_ok=True)

for image in images:
    outfile = outdir / image.with_suffix('.subim').name
    _ = imhead(f'{image}', mode='summary')
    imsubimage(f'{image}', outfile=f'{outfile}', box=box)
    fitsimage = outfile.with_suffix('.subim.fits')
    exportfits(imagename=f'{outfile}', fitsimage=f'{fitsimage}')
