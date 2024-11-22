import sys
sys.path.append('./')
from selfcal_helpers import rank_refants
sys.path.append('/home/users/folguin/alma/analysis_scripts/')
import analysisUtils as au

field = 'G336.018-0.827'
flagchannelsc8 = "0:0~20;130~158;256~274;416~465;545~569;748~774;1355~1368;1610~1632;1835~1850;2361~2386;2399~2411;3830~3839,1:0~10;59~87;389~416;511~538;574~592;1669~1708;1733~1747;1898~1910;1964~1987;2389~2407;2466~2477;2543~2571;2839~2870;3591~3625;3830~3839,2:0~10;28~44;229~266;451~521;523~548;603~646;756~791;968~997;1055~1068;1656~1681;2158~2172;2861~2886;2983~3011;3541~3563;3584~3615;3830~3839,3:0~10;300~371;378~423;440~475;487~500;519~557;605~672;690~700;731~764;863~885;1001~1033;1036~1046;1574~1596;1675~1705;1935~1969;2241~2280;2364~2416;2437~2455;2501~2524;2537~2568;2765~2789;2852~2865;3067~3086;3523~3558;3830~3839"
#flagchannels = "0:,1:,2:,3:"
widths_avgc8 = [120,120,120,120]
widths_avgc9 = [64, 64, 64, 64]
imsize = 10800
cell = '0.004arcsec'
spw = '0,1,2,3'

# Vis
visc8 = f'../uvdata/{field}.c8.shifted.fixvis.ms'
visc9 = f'../uvdata/{field}.c9.ms'
vis_cont_avg_c9 = f'../uvdata/{field}.cont_avg.c9.ms'
vis_cont_avg_c8 = f'../uvdata/{field}.cont_avg.c8.ms'
vis_cont_avg = f'../uvdata/{field}.cont_avg.ms'

# C9
flagmanager(vis=visc9, mode='save', versionname='before_cont_flags')
initweights(vis=visc9, wtmode='weight', dowtsp=True)
#flagdata(vis=vis, mode='manual', spw=flagchannels, flagbackup=False)
split(vis=visc9, spw=spw, outputvis=vis_cont_avg_c9, width=widths_avgc9,
      datacolumn='data')
flagmanager(vis=visc9, mode='restore', versionname='before_cont_flags')

# C8
flagmanager(vis=visc8, mode='save', versionname='before_cont_flags')
initweights(vis=visc8, wtmode='weight', dowtsp=True)
flagdata(vis=visc8, mode='manual', spw=flagchannelsc8, flagbackup=False)
split(vis=visc8, spw=spw, outputvis=vis_cont_avg_c8, width=widths_avgc8,
      datacolumn='data')
flagmanager(vis=visc8, mode='restore', versionname='before_cont_flags')

# Concat
concat(vis=[vis_cont_avg_c9, vis_cont_avg_c8], concatvis=vis_cont_avg)
spw = '0,1,2,3,4,5,6,7'

# Refant from combined
_ = listobs(vis_cont_avg)
rank_refants(vis_cont_avg)
au.commonAntennas(vis_cont_avg)
refant = 'PM04,DV13,DA52,DV08,DV03'

# Correct coordinate
model = f'../../manual_selfcal/G336.018-00.827.4.cont.model'
# Phase cal
# Step 0 --> dirty
threshold = '4mJy'
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.0.cont',
       startmodel=model,
       imsize=[imsize,imsize],
       spw=spw,
       cell=cell,
       specmode='mfs', 
       outframe='LSRK', 
       gridder='standard', 
       deconvolver='hogbom', 
       interactive=False,
       weighting='briggs', 
       robust=0.5, 
       savemodel='modelcolumn',
       niter=0,
       threshold=threshold)
caltable = f'../manual_selfcal/{field}.coord_correction.cal'
solint = 'inf'
rmtables(caltable)
_ = gaincal(vis=vis_cont_avg,
            caltable=caltable,
            field=field,
            gaintype='T',
            refant=refant,
            combine='scan',
            minsnr=1.5,
            calmode='p',
            solint=solint)
applycal(vis=vis_cont_avg,
         gaintable=[caltable],
         interp="linear")

# Step 1
os.system(f'rm -rf ../manual_selfcal/{field}.1.cont.*')
threshold = '0.15mJy'
flagmanager(vis=vis_cont_avg, mode='save', versionname='before_phase_cal_1')
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.1.cont',
       imsize=[imsize, imsize],
       spw=spw,
       cell=cell,
       specmode='mfs',
       outframe='LSRK',
       gridder='standard',
       deconvolver='hogbom',
       interactive=True,
       weighting='briggs',
       robust=0.5, 
       savemodel='modelcolumn',
       niter=10000,
       threshold=threshold) 
caltable = f'../manual_selfcal/{field}.1.phase.cal'
solint = 'inf'
rmtables(caltable)
_ = gaincal(vis=vis_cont_avg,
            caltable=caltable,
            field=field,
            gaintype='T',
            refant=refant,
            minsnr=1.5,
            calmode='p',
            solint=solint,
            combine='spw')
applycal(vis=vis_cont_avg,
         gaintable=[caltable],
         spwmap=['0','0','0','0','4','4','4','4'],
         interp='linear')

# Step 2
os.system(f'rm -rf ../manual_selfcal/{field}.2.cont.*')
threshold = '0.08mJy'
flagmanager(vis=vis_cont_avg, mode='save', versionname='before_phase_cal_2')
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.2.cont',
       imsize=[imsize, imsize],
       spw=spw,
       cell=cell,
       specmode='mfs', 
       outframe='LSRK', 
       gridder='standard', 
       deconvolver='hogbom', 
       interactive=True,
       weighting='briggs', 
       robust=0.5, 
       savemodel='modelcolumn',
       niter=10000,
       threshold=threshold) 
caltable = f'../manual_selfcal/{field}.2.phase.cal'
solint = '45s'
rmtables(caltable)
_ = gaincal(vis=vis_cont_avg,
            caltable=caltable,
            field=field,
            gaintype='T',
            refant=refant,
            minsnr=1.2,
            calmode='p',
            solint=solint,
            combine='spw')
applycal(vis=vis_cont_avg,
         gaintable=[caltable],
         spwmap=['0','0','0','0','4','4','4','4'],
         interp='linear')

# Step 3
os.system(f'rm -rf ../manual_selfcal/{field}.3.cont.*')
threshold = '0.05mJy'
flagmanager(vis=vis_cont_avg, mode='save', versionname='before_phase_cal_3')
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.3.cont',
       imsize=[imsize, imsize],
       spw=spw,
       cell=cell,
       specmode='mfs', 
       outframe='LSRK', 
       gridder='standard', 
       deconvolver='hogbom', 
       interactive=True,
       weighting='briggs', 
       robust=0.5, 
       savemodel='modelcolumn',
       niter=10000,
       threshold=threshold) 
caltable = f'../manual_selfcal/{field}.3.phase.cal'
solint = '20s'
rmtables(caltable)
_ = gaincal(vis=vis_cont_avg,
            caltable=caltable,
            field=field,
            gaintype='T',
            refant=refant,
            minsnr=1.5,
            calmode='p',
            solint=solint,
            combine='spw')
applycal(vis=vis_cont_avg,
         gaintable=[caltable],
         spwmap=['0','0','0','0','4','4','4','4'],
         interp='linear')

# Step 4
os.system(f'rm -rf ../manual_selfcal/{field}.4.cont.*')
threshold = '0.025mJy'
flagmanager(vis=vis_cont_avg, mode='save', versionname='before_amp_cal_1')
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.4.cont',
       imsize=[imsize, imsize],
       spw=spw,
       cell=cell,
       specmode='mfs', 
       outframe='LSRK', 
       gridder='standard', 
       deconvolver='hogbom', 
       interactive=True,
       weighting='briggs', 
       robust=0.5, 
       savemodel='modelcolumn',
       niter=10000,
       threshold=threshold) 

# Step 4 Amp cal
caltable = f'../manual_selfcal/{field}.4.amp.cal'
solint = '45s'
rmtables(caltable)
_ = gaincal(vis=vis_cont_avg,
            caltable=caltable,
            field=field,
            gaintype='T',
            minsnr=1.5,
            refant=refant,
            calmode='ap',
            solnorm=True,
            solint=solint,
            combine='spw')
applycal(vis=vis_cont_avg,
         gaintable=[caltable],
         spwmap=['0','0','0','0','4','4','4','4'],
         interp='linear')

# Final cleaning
os.system(f'rm -rf ../manual_selfcal/{field}.final.cont.*')
threshold = '0.025mJy'
flagmanager(vis=vis_cont_avg, mode='save', versionname='after_amp_cal')
tclean(vis=vis_cont_avg,
       field=field,
       imagename=f'../manual_selfcal/{field}.final.cont',
       mask=f'../manual_selfcal/{field}.4.cont.mask',
       imsize=imsize,
       spw=spw,
       cell=cell,
       specmode='mfs', 
       outframe='LSRK', 
       gridder='standard', 
       deconvolver='hogbom', 
       #interactive=True,
       weighting='briggs', 
       robust=0.5, 
       niter=100000,
       threshold=threshold)
