import numpy as np
from astropy.io import fits
from astropy.wcs.utils import skycoord_to_pixel
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from K2fov import fields

def maketpf(filelist,ra,dec,tpfsize):

 files = np.loadtxt(filelist,dtype=str)
# Make size odd so that target is centered:
 if (tpfsize % 2 == 0): tpfsize += 1 

# Get detector coordinates to cut out around
 midnum = int(len(files)/2)-1
 instr = fits.open(files[midnum],mode='readonly',memmap=True)
 campaign = instr[0].header['CAMPAIGN']     
 instr.close()

 outfile = 'ktwo-'+str(np.round(ra,4))+'-'+str(np.round(dec,4))+'-c'+str(campaign)+'_lpd-targ.fits'

 fovobj = fields.getKeplerFov(campaign)
 ch, detcol, detrow = fovobj.getChannelColRow(ra, dec)
 detcol = np.round(detcol)
 detrow = np.round(detrow)

 skyposition = SkyCoord(ra, dec, unit=('deg','deg'), frame='icrs')
 wcs = WCS(files[midnum])
 pixelpos = skycoord_to_pixel(skyposition, wcs=wcs)
 column = int(np.round(pixelpos[0]))
 row = int(np.round(pixelpos[1]))

# For now, use a template to write the new TPF:
 template = '/Volumes/Work/Field_5/M67/tpf_template.fits'
 tpf = fits.open(template,mode='readonly',memmap=True) 

 timearray = np.empty(len(files),dtype=float)
 timecorrarray = np.empty(len(files),dtype=float)
 cadencearray = np.empty(len(files),dtype=float)
 qualityarray = np.empty(len(files),dtype=float)
 pixels = np.empty([len(files),tpfsize,tpfsize])

 for i in np.arange(len(files)):
   instr = fits.open(files[i],mode='readonly',memmap=True)
   allpixels = instr[0].data[:]
   btlx = instr[0].header['CRVAL1P']
   btly = instr[0].header['CRVAL2P']
   dimx = instr[0].header['NAXIS1']
   dimy = instr[0].header['NAXIS2']
   campaign = instr[0].header['CAMPAIGN']

   timearray[i] = instr[0].header['MIDTIME'] - 2454833.
   timecorrarray[i] = instr[0].header['TIMECORR']
   cadencearray[i] = instr[0].header['CADENCEN']
   qualityarray[i] = instr[0].header['QUALITY']
   pixels[i,:] = allpixels[row-int((tpfsize-1)/2):row+int((tpfsize-1)/2)+1,column-int((tpfsize-1)/2):column+int((tpfsize-1)/2)+1]
   instr.close()

 pixels = np.reshape(pixels,(len(files),tpfsize * tpfsize))

 cards0 = tpf[0].header.cards
 cards1 = tpf[1].header.cards
 cards2 = tpf[2].header.cards

 # construct output primary extension
 hdu0 = fits.PrimaryHDU()
 for i in range(len(cards0)):
        try:
            if cards0[i].keyword not in hdu0.header.keys():
                hdu0.header[cards0[i].keyword] = (cards0[i].value, cards0[i].comment)
            else:
                hdu0.header.cards[cards0[i].keyword].comment = cards0[i].comment
        except:
            pass

# Change a few header keywords:
 hdu0.header['OBJECT']=''
 hdu0.header['KEPLERID']=''
 hdu0.header['CHANNEL']=''
 hdu0.header['CHANNEL']= int(ch)
 hdu0.header['OUTPUT']=''
 hdu0.header['RA_OBJ'] =ra
 hdu0.header['DEC_OBJ'] = dec

 outstr = fits.HDUList(hdu0)


# construct output light curve extension
 coldim = '(' + str(tpfsize) + ',' + str(tpfsize) + ')'
 eformat = str(tpfsize*tpfsize) + 'E'
 jformat = str(tpfsize*tpfsize) + 'J'
 kformat = str(tpfsize*tpfsize) + 'K'
 col1 =  fits.Column(name='TIME', format='D', unit='BJD - 2454833',
                         array=timearray)
 col2 =  fits.Column(name='TIMECORR', format='E', unit='d',
                          array=timecorrarray)
 col3 =  fits.Column(name='CADENCENO', format='J', array=cadencearray)
 col5 =  fits.Column(name='FLUX', format=eformat, unit='e-/s', dim=coldim,
                          array=pixels)
 col10 = fits.Column(name='QUALITY', format='J', array=qualityarray)
 cols =  fits.ColDefs([col1, col2, col3, col5, col10])
 hdu1 =  fits.BinTableHDU.from_columns(cols)

 for i in range(len(cards1)):
        try:
            if cards1[i].keyword not in hdu1.header.keys():
                hdu1.header[cards1[i].keyword] = (cards1[i].value,
                                                  cards1[i].comment)
            else:
                hdu1.header.cards[cards1[i].keyword].comment = cards1[i].comment
        except:
            pass


 hdu1.header['1CRV4P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRV4P'] = detcol - (tpfsize-1)/2 
 hdu1.header['1CRPX4'] = (tpfsize + 1)/2
 hdu1.header['2CRPX4'] = (tpfsize + 1)/2
 hdu1.header['1CRV5P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRV5P'] = detcol - (tpfsize-1)/2
 hdu1.header['1CRPX5'] = (tpfsize + 1)/2
 hdu1.header['2CRPX5'] = (tpfsize + 1)/2
 hdu1.header['1CRV6P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRV6P'] = detcol - (tpfsize-1)/2
 hdu1.header['1CRPX6'] = (tpfsize + 1)/2
 hdu1.header['2CRPX6'] = (tpfsize + 1)/2
 hdu1.header['1CRV7P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRV7P'] = detcol - (tpfsize-1)/2
 hdu1.header['1CRPX7'] = (tpfsize + 1)/2
 hdu1.header['2CRPX7'] = (tpfsize + 1)/2
 hdu1.header['1CRV8P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRP8P'] = detcol - (tpfsize-1)/2
 hdu1.header['1CRPX8'] = (tpfsize + 1)/2
 hdu1.header['2CRPX8'] = (tpfsize + 1)/2
 hdu1.header['1CRV9P'] = detrow - (tpfsize-1)/2
 hdu1.header['2CRV9P'] = detcol - (tpfsize-1)/2
 hdu1.header['1CRPX9'] = (tpfsize + 1)/2
 hdu1.header['2CRPX9'] = (tpfsize + 1)/2

 hdu1.header['NAXIS2'] = len(files)


 outstr.append(hdu1)

# Write output file
 outstr.writeto(outfile,checksum=True)

