"""Sub-images for online storage."""
from pathlib import Path
from casatasks import imsubimage, exportfits, imhead

outdir = Path('../data/cubes/subimages')
indir = Path('../data/cubes/')
images = [
    indir / 'G336.018-0.827_spw0_CH3OH_robust2_multiscale.image.fits',
    indir / 'G336.018-0.827_spw2_SO_robust2_multiscale.image.fits',
    indir / 'G336.018-0.827_spw3_SO_robust2_multiscale.image.fits']
box = '2700,2700,8100,8100'
outdir.mkdir(parents=True, exist_ok=True)

for image in images:
    outfile = outdir / image.with_suffix('.subim.image').name
    _ = imhead(f'{image}', mode='summary')
    imsubimage(f'{image}', outfile=f'{outfile}', box=box)
    fitsimage = outfile.with_suffix('.image.fits')
    exportfits(imagename=f'{outfile}', fitsimage=f'{fitsimage}')
