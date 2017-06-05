import os
import sys
from time import time
import math
import subprocess
from glob import glob
import pandas as pd
from collections import OrderedDict as OD
import easyaccess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.io import fits
from astropy.table import Table
import fitsio
import psycopg2
import HTML

def prep_environ(rootdir,indir,outdir,season,setupfile,version_hostmatch,db,schema):
    os.environ['ROOTDIR']=rootdir
    os.environ['ROOTDIR2']=outdir
    os.environ['INDIR']=indir
    os.environ['EXPDIR']=os.path.join(rootdir,'exp')
    os.environ['SEASON']=season
    os.environ['SETUPFILE']=setupfile
    os.environ['TOPDIR_HOSTMATCH'] = os.path.join(rootdir,'hostmatch/'+version_hostmatch)
    os.environ['OUTDIR_HOSTMATCH'] = os.path.join(outdir,'hostmatch')
    os.environ['DB'] = db
    os.environ['SCHEMA'] = schema

def masterlist(filename,blacklist_file,seqid,propid,expnums=None,a_blacklist=None):
    indir = os.environ.get('INDIR')
    outdir = os.environ.get('ROOTDIR2')
    outdir = os.path.join(outdir,'masterlist')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    
    filename = os.path.join(outdir,filename)

    if os.path.isfile(os.path.join(indir,blacklist_file)):
        blacklist = list(np.genfromtxt(blacklist_file,usecols=(0),unpack=True))
    else:
        blacklist = []

    if a_blacklist:
        blacklist = np.concatenate((blacklist,a_blacklist),axis=0)

    blacklist = [int(x) for x in blacklist]

    if expnums:
        query_exp = """select id as expnum, ra, declination as dec, filter, exptime, airmass, seeing, qc_teff, seqnum, program, object as hex, EXTRACT(EPOCH FROM date - '1858-11-17T00:00:00Z')/(24*60*60) as mjd, TO_CHAR(date - '12 hours'::INTERVAL, 'YYYYMMDD') AS nite 
from exposure 
where propid='2016B-0124' and ra is not null 
and id IN """+str(tuple(expnums))+""" order by id"""
        query_count = """select * from (
WITH objnights AS (
SELECT obstac.nightmjd(date), object, ra, declination
FROM exposure.exposure
WHERE delivered
      AND propid='2016B-0124'
      AND seqid="""+seqid+"""
      AND id IN """+str(tuple(expnums))+"""
GROUP BY obstac.nightmjd(date), object,ra,declination
)
SELECT COUNT(*), ra, declination as dec, object as hex
FROM objnights
GROUP BY object,ra,declination
) as foo order by ra""" 

    else:
        query_exp = """select id as expnum, ra, declination as dec, filter, exptime, airmass, seeing, qc_teff, seqnum, program, object as hex, EXTRACT(EPOCH FROM date - '1858-11-17T00:00:00Z')/(24*60*60) as mjd, TO_CHAR(date - '12 hours'::INTERVAL, 'YYYYMMDD') AS nite 
from exposure 
where propid='2016B-0124' and ra is not null 
order by id"""
        query_count = """select * from (
WITH objnights AS (
SELECT obstac.nightmjd(date), object, ra, declination
FROM exposure.exposure
WHERE delivered
      AND propid='2016B-0124'
      AND seqid="""+seqid+"""
GROUP BY obstac.nightmjd(date), object, ra,declination
)
SELECT COUNT(*), ra, declination as dec, object as hex
FROM objnights
GROUP BY object,ra,declination
) as foo order by ra"""

    conn =  psycopg2.connect(database='decam_prd',
                               user='decam_reader',
                               host='des20.fnal.gov',
                               #password='THEPASSWORD',
                               port=5443) 

    print query_exp
    #print
    #print query_count

    expdf = pd.read_sql(query_exp,conn)

    #ctdf = pd.read_sql(query_count,conn)

    conn.close()

    #print
    #print list(expdf)

    expdf = expdf.loc[~expdf['expnum'].isin(blacklist)]

    expdf = expdf.sort_values(by=['ra','mjd'])

    expdf['dup'] = expdf.duplicated(subset=['ra','nite'])

    epoch = []

    noDupes = []
    [noDupes.append(i) for i in expdf['ra'] if not noDupes.count(i)]
    for x in noDupes:
        ep = 0
        for y in range(len(expdf['ra'])):
            if x==expdf['ra'][y] and expdf['dup'][y]==True:
                ep = ep
                epoch.append(ep)
            elif x==expdf['ra'][y]:
                ep = ep + 1
                epoch.append(ep)

    ### strip hex string to just ra/dec term
    striphex = []
    for ihex in expdf['hex']:
        a = ihex
        a = a.split('hex')
        a = a[1].split('tiling')
        new = a[0].strip()
        striphex.append(new)
    
    niteform = lambda x: int(x)
    expdf['nite'] = expdf['nite'].map(niteform)

    mjdform = lambda x: round(x,3)
    expdf['mjd'] = expdf['mjd'].map(mjdform)

    expdf['epoch'] = epoch

    expdf['striphex'] = striphex
    
    tbhdu1 = fits.BinTableHDU.from_columns(
        [fits.Column(name='fullhex', format='A69', array=expdf['hex']),
         fits.Column(name='hex', format='A8', array=striphex),
         fits.Column(name='epoch', format='K', array=expdf['epoch']),
         fits.Column(name='expnum', format='K', array=expdf['expnum']),
         fits.Column(name='RA', format='E', array=expdf['ra']),
         fits.Column(name='DEC', format='E', array=expdf['dec']),
         fits.Column(name='nite', format='K', array=expdf['nite']),
         fits.Column(name='mjd', format='E', array=expdf['mjd']),
         fits.Column(name='t_eff', format='E', array=expdf['qc_teff']),
         ])

    tbhdu1.writeto(filename,clobber=True)
    
    f=open(os.path.join(outdir,blacklist_file),'w')
    f.write(str(sorted(set(blacklist))))
    f.close()
    
    return expdf[['expnum','nite','filter']],filename

def checkoutputs(expdf,logfile,ccdfile,goodchecked,steplist):
    expnums = expdf['expnum'].tolist()
    nites = expdf['nite'].tolist()
    bands = expdf['filter'].tolist()

    season = os.environ.get('SEASON')
    outdir = os.path.join(os.environ.get('ROOTDIR2'),'checkoutputs')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    steplist = os.path.join(os.environ.get('INDIR'),steplist)
    f = open('steplist.txt','r')
    stepnames = f.readlines()
    f.close()

    stepnames = map(lambda x: x.strip(), stepnames)

    goodchecked = os.path.join(outdir,goodchecked)
    if os.path.isfile(goodchecked):
        f = open(goodchecked,'r')
        good = f.readlines()
        f.close()
        good = map(lambda x: int(x.strip()), good)        
    else:
        good = []

    expdir = os.environ.get('EXPDIR')
    logname = os.path.join(outdir,logfile)
    lf = open(logname,'w+')
    lf.write('EXPOSURES PROVIDED: '),lf.write(','.join(map(str,sorted(expnums))))
    lf.write('\n\n')
    d = OD()
    d['expnum'] = []
    chips = range(1,63)
    steps = range(1,29)
    chips.remove(2),chips.remove(31),chips.remove(61)
    for ch in chips:
        ch = '%02d' % ch
        d[ch] = []
    for e in expnums:
        nite = nites[expnums.index(e)]
        nite = str(nite)
        band = bands[expnums.index(e)]
        e = str(e)
        end = nite+'/'+e+'/'+'dp'+season
        p = os.path.join(expdir,end)
        if int(e) not in good:
            if os.path.isdir(p):
                print str(expnums.index(int(e)))+'/'+str(len(expnums))+' - '+p
                d['expnum'].append(e)
            else:
                print str(expnums.index(int(e)))+'/'+str(len(expnums))+' - '+p, 'does not exist. Check diffimg outputs.'
                continue
            for c in chips:
                c = '%02d' % c
                p2 = os.path.join(p,band+'_'+c)
                if os.path.isdir(p2):
                    for r in steps:
                        fail = stepnames[r-1]+'.FAIL'
                        gpfail = os.path.join(p2,fail)
### The current assumption is that .FAIL files are cleared out when a CCD is reprocessed. 
### If this is not true, uncomment the 3 lines below and tab the append and break lines. 
### In that event, one must also consider how to deal with a RUN28 failure.
                        if os.path.isfile(gpfail):
                            #log = stepnames[r]+'.LOG'
                            #plog = os.path.join(gp,log)
                            #if os.path.isfile(plog):
                            d[c].append(int(r))
                            break
                    else:
                        d[c].append(0)
        else:
            print str(expnums.index(int(e)))+'/'+str(len(expnums))+' - '+p
            d['expnum'].append(e)
            for c in chips:
                c = '%02d' % c
                d[c].append(0)

    lf.write('EXPOSURES CHECKED: '),lf.write(','.join(map(str,sorted(d['expnum']))))
    lf.write('\n\n')

    nonex,yesex = [],[]
    for x in sorted(expnums):
        x=str(x)
        if x not in d['expnum']:
            nonex.append(x)
        else:
            yesex.append(x)
    lf.write('EXPOSURES NOT FOUND: ')
    if len(nonex)==0:
        lf.write('none')
    else:
        lf.write(','.join(map(str,nonex)))
    lf.write('\n\n')
    
    df1 = pd.DataFrame(d)
    df = df1.set_index('expnum')

    ccddf = df.copy()

    listgood = df.loc[df.sum(axis=1) == 0].index
    listgood = listgood.tolist()
    listgood = map(lambda x: int(x), listgood)
    np.savetxt(goodchecked,sorted(listgood),fmt='%d')

    df['successes']=(df==0).astype(int).sum(axis=1)
    df['fraction'] = ""
    for exp in list(df.index.values):
        frac = float(df.get_value(exp,'successes'))/59.
        frac = round(frac,3)
        df.set_value(exp,'fraction',frac)
    
    lf.write('# OF CCDS PROCESSED: ')
    ccdsum = sum(df['successes'])
    lf.write(str(ccdsum)),lf.write('\n')
    
    lf.write('TOTAL CCDS CHECKED: ')
    ccdtot = 59*len(df['successes'])
    lf.write(str(ccdtot)),lf.write('\n')
    
    lf.write('CCD SUCCESS RATE: ')
    div = float(ccdsum)/float(ccdtot)
    lf.write(str(round(div,3)))
    lf.close()
    
    df.sort_index(inplace=True)
    df.to_csv(os.path.join(outdir,ccdfile))

    return yesex,nonex,ccddf
            
def forcephoto(ncore=4,numepochs_min=0,writeDB=False):    
    season = os.environ.get('SEASON')
    a = './forcePhoto_master.pl ' 
    a = a + ' -season ' + season 
    a = a + ' -numepochs_min ' + numepochs_min 
    a = a + ' -ncore ' + ncore 
    a = a + ' -noprompt ' 
    if writeDB == True:
        a = a + ' -writeDB ' 
    print a
    a = 'source '+os.getenv('SETUPFILE')+'; '+a
    #subprocess.call(a,shell=True)
 
def truthtable(expnums,filename,truthplus):
    season = os.environ.get('SEASON')
    outdir = os.path.join(os.environ.get('ROOTDIR2'),'truthtable')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    db = os.environ.get('DB')
    schema = os.environ.get('SCHEMA')

    explist=','.join(map(str,expnums))

### Truth table (normal)
    query='select distinct SNFAKE_ID, EXPNUM, CCDNUM, TRUEMAG, TRUEFLUXCNT, FLUXCNT, BAND, NITE, MJD, SEASON from '+ schema +'.SNFAKEIMG where EXPNUM IN ('+explist+') and SEASON='+ season +' order by SNFAKE_ID'
    print query

    filename=os.path.join(outdir,filename)
    connection=easyaccess.connect(db)
    connection.query_and_save(query,filename)

    print

### Truth table plus

    #query='select f.SNFAKE_ID, f.EXPNUM, f.CCDNUM, o.RA, o.DEC, o.MAG, o.FLUX, o.FLUX_ERR, f.TRUEMAG, f.TRUEFLUXCNT, o.FLUX, o.SEXFLAGS, f.BAND, f.NITE, f.MJD, f.SEASON from '+ schema +'.SNFAKEIMG f, '+ schema +'.SNOBS o where f.SNFAKE_ID=o.SNFAKE_ID and f.EXPNUM=o.EXPNUM and f.SEASON='+ season +' and f.SEASON=o.SEASON order by SNFAKE_ID'

    query = 'select SNFAKE_ID, EXPNUM, CCDNUM, RA, DEC, -2.5*log(10,FLUXCNT)+ZERO_POINT as MAG, MAGOBS_ERR as MAGERR, FLUXCNT, TRUEMAG, TRUEFLUXCNT, SNR_DIFFIM as SNR, REJECT, ML_SCORE, BAND, NITE, SEASON from '+ schema +'.SNFAKEMATCH where SEASON='+ season +' order by SNFAKE_ID'

    print query

    plus = connection.query_to_pandas(query)
    connection.query_and_save(query,os.path.join(outdir,truthplus))

    connection.close()

    return plus

def makedatafiles(format,numepochs_min,two_nite_trigger,outfile,outdir,fakeversion=None):
    season = os.environ.get('SEASON')
    datafiles_dir = os.path.join(os.environ.get('ROOTDIR2'),'makedatafiles')
    if not os.path.isdir(datafiles_dir):
        os.mkdir(datafiles_dir)
    a = 'makeDataFiles_fromSNforce' 
    a = a + ' -format ' + format 
    a = a + ' -season ' + season  
    a = a + ' -numepochs_min ' + numepochs_min  
    a = a + ' -outFile_stdout ' + outfile 
    a = a + ' -outDir_data ' + outdir
    if not two_nite_trigger == 'null':
        a = a + ' -2nite_trigger ' + trigger 
    if not fakeversion == None: 
        a = a + ' -fakeVersion ' + fakeversion
    a = 'source '+os.getenv('SETUPFILE')+'; cd '+datafiles_dir+ '; '+ a + '; cd -'
    print a
    subprocess.call(a, shell=True)
    
def combinedatafiles(master,fitsname,datadir):
    season = os.environ.get('SEASON')
    season = str(season)

    mlist = Table.read(master)
    masdf = mlist.to_pandas()

    path = os.path.join(os.environ.get('ROOTDIR2'), 'makedatafiles')
    fitsname = os.path.join(path,fitsname)
    path = os.path.join(path,datadir)
    
    if os.path.isfile(fitsname):
        print 'A combined .fits file for all real candidates already exists in the specified outdir with the specified name:'
        print
        print fitsname
        print
        print 'If you want to recreate the file, either change the combined_fits key under the [GWmakeDataFiles-real] heading in the .ini file, or simply delete the existing one.'
        print
        return fitsname

    listfile = os.path.join(path,datadir+'.LIST')

    ls = open(listfile,'r+')
    dats = ls.readlines()
    ls.close()

    hostlist = []
    c = 0

    MJD,BAND,FIELD,FLUXCAL,FLUXCALERR,PHOTFLAG,PHOTPROB,ZPFLUX,PSF,SKYSIG,\
        SKYSIG_T,GAIN,XPIX,YPIX,NITE,EXPNUM,CCDNUM,OBJID = [],[],[],[],[],[],\
        [],[],[],[],[],[],[],[],[],[],[],[]

    RA,DEC,CAND_ID,DATAFILE,SN_ID = [],[],[],[],[]
    HOSTID,PHOTOZ,PHOTOZERR,SPECZ,SPECZERR,HOSTSEP,HOST_GMAG,HOST_RMAG,HOST_IMAG,\
        HOST_ZMAG = [],[],[],[],[],[],[],[],[],[]

    c=0
    allgood=0
    for d in dats:
        c=c+1
        if c%1000==0:
            print c
            #break                                                                                         
        filename = d.split('\n')[0]
        datfile = os.path.join(path,filename)
        f = open(datfile,'r+')
        lines = f.readlines()
        f.close()

        snid = lines[1].split()[1]
        raval = lines[8].split()[1]
        decval = lines[9].split()[1]
        host_id = lines[15].split()[1]
        photo_z = lines[16].split()[1]
        photo_zerr = lines[16].split()[3]
        spec_z = lines[17].split()[1]
        spec_zerr = lines[17].split()[3]
        host_sep = lines[18].split()[1]
        h_gmag = lines[19].split()[1]
        h_rmag = lines[19].split()[2]
        h_imag = lines[19].split()[3]
        h_zmag = lines[19].split()[4]

        mjd,band,field,fluxcal,fluxcalerr,photflag,photprob,zpflux,psf,skysig,skysig_t,gain,xpix,ypix,nite,expnum,ccdnum,objid = np.genfromtxt(datfile,skip_header=53,usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),unpack=True)

        band = np.genfromtxt(datfile,dtype='string',skip_header=53,usecols=(2,),unpack=True)

        if all(x==12288 for x in photflag):
            allgood=allgood+1

        n = len(mjd)

        ra = np.empty(n)
        ra.fill(raval)
        dec = np.empty(n)
        dec.fill(decval)
        cand = np.empty(n)
        cand.fill(c)
        sn_id = np.empty(n)
        sn_id.fill(snid)
        hostid = np.empty(n)
        hostid.fill(host_id)
        photoz=np.empty(n)
        photoz.fill(photo_z)
        photozerr=np.empty(n)
        photozerr.fill(photo_zerr)
        specz=np.empty(n)
        specz.fill(spec_z)
        speczerr=np.empty(n)
        speczerr.fill(spec_zerr)
        hostsep=np.empty(n)
        hostsep.fill(host_sep)
        hgmag=np.empty(n)
        hgmag.fill(h_gmag)
        hrmag=np.empty(n)
        hrmag.fill(h_rmag)
        himag=np.empty(n)
        himag.fill(h_imag)
        hzmag=np.empty(n)
        hzmag.fill(h_zmag)

        for j in range(n):
            DATAFILE.append(filename)

        for k in range(n):
            RA.append(ra[k])
            DEC.append(dec[k])
            CAND_ID.append(cand[k])
            SN_ID.append(sn_id[k])

            HOSTID.append(hostid[k])
            PHOTOZ.append(photoz[k])
            PHOTOZERR.append(photozerr[k])
            SPECZ.append(specz[k])
            SPECZERR.append(speczerr[k])
            HOSTSEP.append(hostsep[k])
            HOST_GMAG.append(hgmag[k])
            HOST_RMAG.append(hrmag[k])
            HOST_IMAG.append(himag[k])
            HOST_ZMAG.append(hzmag[k])

            MJD.append(mjd[k])
            BAND.append(band[k])
            FIELD.append(field[k])
            FLUXCAL.append(fluxcal[k])
            FLUXCALERR.append(fluxcalerr[k])
            PHOTFLAG.append(photflag[k])
            PHOTPROB.append(photprob[k])
            ZPFLUX.append(zpflux[k])
            PSF.append(psf[k])
            SKYSIG.append(skysig[k])
            SKYSIG_T.append(skysig_t[k])
            GAIN.append(gain[k])
            XPIX.append(xpix[k])
            YPIX.append(ypix[k])
            NITE.append(nite[k])
            EXPNUM.append(expnum[k])
            CCDNUM.append(ccdnum[k])
            OBJID.append(objid[k])

    print
    #print len(RA)
    #print len(PHOTFLAG)

    HEX = []

    for h in EXPNUM:
        HEX.append(masdf['hex'].loc[masdf['expnum']==h].values[0])

    MJD,FIELD,FLUXCAL,FLUXCALERR,PHOTFLAG,PHOTPROB,ZPFLUX,PSF,SKYSIG,SKYSIG_T,\
        GAIN,XPIX,YPIX,NITE,EXPNUM,CCDNUM,OBJID,RA,DEC,CAND_ID,SN_ID = \
        np.asarray(MJD),np.asarray(FIELD),np.asarray(FLUXCAL),np.asarray(FLUXCALERR),\
        np.asarray(PHOTFLAG),np.asarray(PHOTPROB),np.asarray(ZPFLUX),np.asarray(PSF),\
        np.asarray(SKYSIG),np.asarray(SKYSIG_T),np.asarray(GAIN),np.asarray(XPIX),\
        np.asarray(YPIX),np.asarray(NITE),np.asarray(EXPNUM),np.asarray(CCDNUM),\
        np.asarray(OBJID),np.asarray(RA),np.asarray(DEC),np.asarray(CAND_ID),np.asarray(SN_ID)

    HOSTID,PHOTOZ,PHOTOZERR,SPECZ,SPECZERR,HOSTSEP,HOST_GMAG,HOST_RMAG,HOST_IMAG,\
        HOST_ZMAG = np.asarray(HOSTID),np.asarray(PHOTOZ),np.asarray(PHOTOZERR),\
        np.asarray(SPECZ),np.asarray(SPECZERR),np.asarray(HOSTSEP),np.asarray(HOST_GMAG),\
        np.asarray(HOST_RMAG),np.asarray(HOST_IMAG),np.asarray(HOST_ZMAG)

    tbhdu1 = fits.BinTableHDU.from_columns(
        [fits.Column(name='cand_ID', format='K', array=CAND_ID.astype(float)),
         fits.Column(name='SNID', format='K', array=SN_ID.astype(float)),
         fits.Column(name='OBJID', format='K', array=OBJID.astype(float)),
         fits.Column(name='RA', format='E', array=RA.astype(float)),
         fits.Column(name='DEC', format='E', array=DEC.astype(float)),
         fits.Column(name='MJD', format='E', array=MJD.astype(float)),
         fits.Column(name='BAND', format='1A', array=BAND),
         fits.Column(name='EXPNUM', format='K', array=EXPNUM.astype(float)),
         fits.Column(name='CCDNUM', format='K', array=CCDNUM.astype(float)),
         fits.Column(name='NITE', format='K', array=NITE.astype(float)),
         fits.Column(name='HEX', format='8A', array=HEX),
         #fits.Column(name='FIELD', format='K', array=RA.astype(float)),
         fits.Column(name='FLUXCAL', format='E', array=FLUXCAL.astype(float)),
         fits.Column(name='FLUXCALERR', format='E', array=FLUXCALERR.astype(float)),
         fits.Column(name='PHOTFLAG', format='K', array=PHOTFLAG.astype(float)),
         fits.Column(name='PHOTPROB', format='E', array=PHOTPROB.astype(float)),
         fits.Column(name='ZPFLUX', format='E', array=ZPFLUX.astype(float)),
         fits.Column(name='PSF', format='E', array=PSF.astype(float)),
         fits.Column(name='SKYSIG', format='E', array=SKYSIG.astype(float)),
         fits.Column(name='SKYSIG_T', format='E', array=SKYSIG_T.astype(float)),
         fits.Column(name='GAIN', format='E', array=GAIN.astype(float)),
         fits.Column(name='XPIX', format='E', array=XPIX.astype(float)),
         fits.Column(name='YPIX', format='E', array=YPIX.astype(float)),
    
         fits.Column(name='HOSTID', format='K', array=HOSTID.astype(float)),
         fits.Column(name='PHOTOZ', format='E', array=PHOTOZ.astype(float)),
         fits.Column(name='PHOTOZERR', format='E', array=PHOTOZERR.astype(float)),
         fits.Column(name='SPECZ', format='E', array=SPECZ.astype(float)),
         fits.Column(name='SPECZERR', format='E', array=SPECZERR.astype(float)),
         fits.Column(name='HOSTSEP', format='E', unit='arcsec', array=HOSTSEP.astype(float)),
         fits.Column(name='HOST_GMAG', format='E', array=HOST_GMAG.astype(float)),
         fits.Column(name='HOST_RMAG', format='E', array=HOST_RMAG.astype(float)),
         fits.Column(name='HOST_IMAG', format='E', array=HOST_IMAG.astype(float)),
         fits.Column(name='HOST_ZMAG', format='E', array=HOST_ZMAG.astype(float)),
         fits.Column(name='DATAFILE', format='21A', array=DATAFILE)])

    tbhdu1.writeto(fitsname,clobber=True)

    print "number of candidates where all detections had ml_score>0.5 :",allgood
    print
    
    return fitsname

def makeplots(ccddf,master,truthplus,fitsname,expnums,mjdtrigger,ml_score_cut=0.,skip=False):

    season = os.environ.get('SEASON')
    season = str(season)

    rootdir = os.environ.get('ROOTDIR')
    rootdir = os.path.join(rootdir,'exp')

### get data
    if os.path.isfile(master):
        mlist = Table.read(master)
        masdf = mlist.to_pandas()
    else:
        skip = True
        print "No master list found with filename",master+'.'
        print "Plots requiring a master list (SNR, RA/DEC hex maps) will not be created."

    df1 = truthplus

    outdir = os.path.join(os.environ.get('ROOTDIR2'),'plots')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    rtable = Table.read(fitsname)
    rdf1 = rtable.to_pandas()
    
### cut out rejects ###
    df = df1.loc[df1['REJECT'] == 0]
    rdf = rdf1.loc[rdf1['PHOTFLAG'].isin([4096,12288])]

### TEMPLATE FAILURES (temporary section?) ###
    ccdhexes = masdf['fullhex'].loc[masdf['epoch']==4].values
    ccdexp = masdf.loc[(masdf['epoch']==4) & (masdf['fullhex'].isin(ccdhexes))]

    tempfails,rafail,decfail = [],[],[]

    explist = list(ccdexp['expnum'].values)

    for t in sorted(explist):
        if (ccddf.loc[str(t)]==2).any():
            tempfails.append(t)
            rafail.append(masdf['RA'].loc[masdf['expnum']==t].values[0])
            decfail.append(masdf['DEC'].loc[masdf['expnum']==t].values[0])
    
    print len(tempfails)
    print tempfails

    np.savetxt('event3_ccds.txt',np.c_[tempfails,rafail,decfail],fmt=['%d','%.6f','%.6f'],delimiter='\t',header='EXP\tRA\t\tDEC',comments='')

    notemp = masdf.loc[masdf['expnum'].isin(tempfails)]
    yestemp = masdf.loc[(masdf['expnum'].isin(explist)) & (~masdf['expnum'].isin(tempfails))]

    notemp.ix[notemp['RA'] > 180, 'RA'] = notemp.ix[notemp['RA'] > 180, 'RA'] - 360.
    yestemp.ix[yestemp['RA'] > 180, 'RA'] = yestemp.ix[yestemp['RA'] > 180, 'RA'] - 360.

    xmax = max([max(notemp['RA']),max(yestemp['RA'])])+2
    ymax = max([max(notemp['DEC']),max(yestemp['DEC'])])+2
    xmin = min([min(notemp['RA']),min(yestemp['RA'])])-2
    ymin = max([min(notemp['DEC']),min(yestemp['DEC'])])-2

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    plt.scatter(notemp['RA'],notemp['DEC'],marker='H',c='r',s=200,label='failed')
    plt.scatter(yestemp['RA'],yestemp['DEC'],marker='H',c='b',s=200,label='succeeded')
    plt.title('Template failures - GW170104')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.legend()
    #plt.show()
    plt.savefig('templatefailures.png',dpi=200)
    
    plt.clf()
    sys.exit()
    
### EFFICIENCY ###

    bins = np.arange(17,25,1)

    fhist, bin_edges = np.histogram(df['MAG'], bins=bins)
    thist, bin_edges = np.histogram(df1['TRUEMAG'], bins=bins)
    

    plt.xlim(17,25)
    plt.ylim(0,100)
    plt.plot(bins[:-1], fhist*100.0/thist, lw=4)
    plt.scatter(bins[:-1], fhist*100.0/thist, lw=4)
    plt.title('Efficiency')
    plt.xlabel('Mag')
    plt.ylabel('Percent Found')
    plt.savefig(os.path.join(outdir,'efftest_'+season+'.png'))
    plt.clf()

### ML_SCORE HISTOGRAM - FAKES ###

    plt.hist(df1['ML_SCORE'],bins=np.linspace(0.3,1,100))
    plt.title('ML_SCORE OF FAKES')
    plt.xlabel('ml_score')
    plt.ylabel('# of fakes')
    plt.savefig(os.path.join(outdir,'fakemltest_'+season+'.png'))
    plt.clf()

### PULL --> (MAG-TRUEMAG)/MAG_ERR -- FOR FAKES ### 

    magdiff = df['MAG']-df['TRUEMAG']
    pull = magdiff/df['MAGERR']
    
    #bins = np.linspace(int(min(pull)),int(max(pull)+1),100)
    bins = np.linspace(-3,3,100)
    plt.hist(pull,bins=bins)
    plt.xlabel("Magnitude Pull (MAG-TRUEMAG)/MAGERR")
    plt.ylabel("Number of Objects")
    plt.title("Magnitude Pull for Fakes")
    plt.xlim(bins.min(),bins.max())

    plt.savefig(os.path.join(outdir,'pulltest_'+season+'.png'))
    plt.clf()

### NUMBER OF REAL CANDIDATES PER CCD ###

    bins = np.arange(1,64,1)
    
    for e in expnums:
        ccdcand = rdf['CCDNUM'].loc[rdf['EXPNUM'] == e]
        ccdhist, bin_edges = np.histogram(ccdcand, bins=bins)
        
        plt.xlim(0,63)
        plt.ylim(0,max(ccdhist)+1)
        plt.scatter(bins[:-1], ccdhist)
        plt.plot(bins[:-1], ccdhist, label=str(e))
        plt.grid(True)
    
    plt.xlabel('CCD')
    plt.ylabel('# of candidates')
    plt.title('# of Candidates per CCD')
    plt.legend(fontsize='small')
    plt.savefig(os.path.join(outdir,'ccdtest'+season+'.png'))
    plt.clf()

### SNR VS. HEX -- FAKES ###
    if not skip:
        masdf_ord = masdf
        masdf_ord.ix[masdf_ord['RA'] > 180, 'RA'] = masdf_ord.ix[masdf_ord['RA'] > 180, 'RA'] - 360.
        masdf_ord = masdf_ord.sort_values(by='RA')
        SNR = OD()
        maxepoch = max(masdf_ord['epoch'])
        for hx in masdf_ord['hex'].unique():
            SNR[hx]=np.array([-5.]*maxepoch)
            epexpdf = masdf_ord[['hex','epoch','expnum']].loc[masdf_ord['hex']==hx]
            for ep in epexpdf['epoch'].unique():
                snrexp = epexpdf['expnum'].loc[epexpdf['epoch']==ep].values[0]
                snrs = df['SNR'].loc[df['EXPNUM']==snrexp].values
                if len(snrs)==0:
                    SNR[hx][ep-1] = 0
                else:
                    meansnr = np.mean(snrs)
                    SNR[hx][ep-1] = meansnr

        SNRdf = pd.DataFrame.from_dict(SNR,orient='index')
        SNRdf = SNRdf.reset_index()
        cols = ['hex']

        epochs = np.array(range(maxepoch))+1
        for i in epochs:
            cols.append(str(i))
        SNRdf.columns = cols

        inds = np.array(SNRdf.index.values)+1

        plt.figure(figsize=(16,9))
        
        ax = plt.axes()
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)

        major = MultipleLocator(5)
        majForm = FormatStrFormatter('%d')
        minor = MultipleLocator(1)
        ax.xaxis.set_major_locator(major)
        ax.xaxis.set_major_formatter(majForm)
        ax.xaxis.set_minor_locator(minor)

        for e in epochs:
            plt.plot(inds, SNRdf[str(e)], '-o', markeredgewidth=0, label='epoch '+str(e), antialiased=True)
        
        plt.legend(loc='upper left',fontsize='small')
        plt.axis([0, max(inds)+1, -10, max(SNRdf.max(numeric_only=True))+5])
        plt.title("Average SNR of MAG=20 fakes by hex index and epoch - GW170104")

        plt.xlabel('Hex Index')
        plt.ylabel('Mean SNR')

        savefile = 'SNR_test.png'
        plt.savefig(os.path.join(outdir,savefile),dpi=400)

        plt.clf()

    sys.exit()

### RA/DEC MAPS ###
    
    radecdf = rdf
    if abs(max(radecdf['RA'])-min(radecdf['RA']))>180:
        for ira in range(len(radecdf['RA'])):
            if radecdf['RA'][ira]>180:
                radecdf['RA'][ira] = radecdf['RA'][ira]-360

    radecdf = radecdf.drop_duplicates('cand_ID')

    radecdf = radecdf.loc[radecdf['PHOTPROB'] > ml_score_cut]

    #plt.hist2d(radecdf['RA'],radecdf['DEC'],50)
    
    mapdir = os.path.join(outdir,'maps')
    if not os.path.isdir(mapdir):
        os.mkdir(mapdir)    

    hexex = []

    ### this loop gets the full set of first epoch exposures of each hex.
    ### if there are two (or more), it chooses the one with the best t_eff.
    for h in masdf['fullhex'].unique():
        exepteff = masdf[['expnum','epoch','t_eff']].loc[masdf['fullhex'] == h]
        cut = exepteff[['expnum','epoch','t_eff']].loc[exepteff['epoch']==1]
        if len(cut)>1:
            cut = cut.loc[cut['t_eff'] == cut['t_eff'].ix[cut['t_eff'].idxmax()]]
        hexex.append(cut['expnum'].values[0])

    radecdf = radecdf.loc[radecdf['EXPNUM'].isin(hexex)]

    ### overall map
    plt.scatter(radecdf['RA'],radecdf['DEC'],c=radecdf['PHOTPROB'],edgecolor='',s=5)
    plt.xlim(min(radecdf['RA'])-0.2,max(radecdf['RA'])+0.2)
    plt.ylim(min(radecdf['DEC'])-0.2,max(radecdf['DEC'])+0.2)
    plt.clim(0,1)
    plt.colorbar().set_label('ml_score')
    plt.title('Candidate Sky Map')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.savefig(os.path.join(outdir,'fullmap_'+season+'.png'))
    plt.clf()

    ### individual hex maps
    if not skip:
        for e in hexex:
            print e
            out = ''
            out = os.path.join(rootdir,str(masdf['nite'].loc[masdf['expnum']==e].values[0]))
            out = os.path.join(out,str(e))
            out = os.path.join(out,str(e)+'.out')

            odf = pd.read_table(out,delim_whitespace=True,header=None,names=['expnum','band','ccd','ra1','dec1','ra2','dec2','ra3','dec3','ra4','dec4'])

            odf = odf.drop_duplicates()
            odf = odf.reset_index(drop=True)

            odf.ix[odf.ra1 > 270., 'ra1'] = odf.ix[odf.ra1 > 270., 'ra1'] - 360.
            odf.ix[odf.ra2 > 270., 'ra2'] = odf.ix[odf.ra2 > 270., 'ra2'] - 360.
            odf.ix[odf.ra3 > 270., 'ra3'] = odf.ix[odf.ra3 > 270., 'ra3'] - 360.
            odf.ix[odf.ra4 > 270., 'ra4'] = odf.ix[odf.ra4 > 270., 'ra4'] - 360.

            ras = np.concatenate((odf['ra1'].tolist(),odf['ra2'].tolist(),odf['ra3'].tolist(),odf['ra4'].tolist()),axis=0)

            decs = np.concatenate((odf['dec1'].tolist(),odf['dec2'].tolist(),odf['dec3'].tolist(),odf['dec4'].tolist()),axis=0)

            for i in range(len(odf)):
                ra = odf.ix[i,['ra1','ra2','ra3','ra4']]
                dec = odf.ix[i,['dec1','dec2','dec3','dec4']]
                chip = str(odf.ix[i,'ccd'])
                midra = (max(ra)+min(ra))/2.
                middec = (max(dec)+min(dec))/2.
                middle = tuple([midra,middec])
                cs = zip(ra,dec)
                cent=(sum([c[0] for c in cs])/len(cs),sum([c[1] for c in cs])/len(cs))
                cs.sort(key=lambda c: math.atan2(c[1]-cent[1],c[0]-cent[0]))
                cs.append(cs[0])
                plt.plot([c[0] for c in cs],[c[1] for c in cs],ls=':',lw=0.5,c='k')
                plt.annotate(chip, xy=middle, ha='center',va='center',family='sans-serif',fontsize=12,alpha=0.3)

            plt.xlim(min(ras)-0.2,max(ras)+0.2)
            plt.ylim(min(decs)-0.2,max(decs)+0.2)

            pltdf = radecdf.loc[radecdf['EXPNUM'] == e]
            plt.scatter(pltdf['RA'],pltdf['DEC'],c=pltdf['PHOTPROB'],edgecolor='')
            plt.clim(0,1)
            plt.colorbar().set_label('ml_score')
            hexname = masdf['hex'].loc[masdf['expnum'] == e].values[0]
            plt.title('Candidate Sky Map: Hex '+str(hexname)+' (Exposure '+str(e)+')')
            plt.xlabel('RA')
            plt.ylabel('DEC')
            plt.savefig(os.path.join(mapdir,'map_'+str(hexname)+'_'+str(e)+'.png'),dpi=200)
            plt.clf()

    lcdir = os.path.join(outdir,'lightcurves')

    band = 'i'

    numsnid = len(rdf['SNID'].unique())
    ctsnid = 0
    noct,yesct = 0,0

    rdf['cutflag'] = np.zeros(len(rdf))

    for sn in rdf['SNID'].unique():
        #ctsnid += 1
        
        new = rdf[['MJD','FLUXCAL','FLUXCALERR','PHOTPROB']].loc[rdf['SNID']==sn]
        
        if all(i < ml_score_cut for i in new['PHOTPROB']):
            noct+=1
            continue
        
        yesct+=1
        
        if (noct+yesct) % 100 == 0:
            print str(noct+yesct)+'/'+str(numsnid)

        rdf.loc[rdf['SNID']==sn,'cutflag'] = 1
        
        continue

        mjd = np.array(new['MJD'].tolist())
        flux = np.array(new['FLUXCAL'].tolist())
        fluxerr = np.array(new['FLUXCALERR'].tolist())
        ml_score = np.array(new['PHOTPROB'].tolist())

        plt.errorbar(mjd-mjdtrigger,flux,yerr=fluxerr,fmt='none',ecolor='k',zorder=0)
        plt.scatter(mjd-mjdtrigger,flux,c=ml_score,edgecolor='',s=40,zorder=1)
        plt.clim(0,1)
        plt.title('SNID '+str(sn)+' ('+band+')')
        plt.colorbar().set_label('ml_score')
        plt.xlabel('MJD - MJD(TRIGGER)')
        plt.ylabel('FLUX')
        plt.savefig(os.path.join(lcdir,'SNID'+str(sn)+'.png'))
        plt.clf()
    return rdf

def createhtml(fitsname,realdf,master):
    rootdir = os.environ.get('ROOTDIR')
    expdir = os.path.join(rootdir,'exp')
    season = os.environ.get('SEASON')
    skip = False

    if os.path.isfile(master):
        mlist = Table.read(master)
        masdf = mlist.to_pandas()
    else:
        skip = True
        print "No master list found with filename",master+'.'
        print "This step will run more slowly because it will require the use of glob."

    rdf = realdf.reset_index(drop=True)

    ### GET STAMPS ###
    lenr = len(rdf)
    srcharray = ['' for x in range(lenr)]
    temparray = ['' for x in range(lenr)]
    diffarray = ['' for x in range(lenr)]
    aaa = 0
    aaalen = len(rdf['EXPNUM'].unique())
    #time1 = time()
    for e in sorted(rdf['EXPNUM'].unique()):
        aaa += 1
        bb = 0
        print str(aaa)+'/'+str(aaalen)+' - '+str(e)
        edf = rdf[['EXPNUM','NITE','CCDNUM','BAND','OBJID','HEX']].loc[rdf['EXPNUM'] == e]
        bblen = len(edf['CCDNUM'].unique())
        #time2 = time()
        for c in sorted(edf['CCDNUM'].unique()):
            bb += 1
            #print '    '+str(bb)+'/'+str(bblen)+' - '+str(c)
            cdf = edf.loc[edf['CCDNUM'] == c]
            nite = str(cdf['NITE'].values[0])
            hhex = str(cdf['HEX'].values[0])
            exp = str(e)
            dp = 'dp'+str(season)
            band = str(cdf['BAND'].values[0])
            ccd = '%02d' % c
            stampname = 'stamps_'+nite+'_'+hhex+'_'+band+'_'+ccd
            stampstar = os.path.join(expdir,nite+'/'+exp+'/'+dp+'/'+band+'_'+ccd+'/'+stampname)
            #time3 = time()
            #gstamp = glob(stampstar)
            #time4 = time()
            #if len(gstamp)==1:
            if os.path.isdir(stampstar):
                stampdir = stampstar
                for i in list(cdf.index.values):
                    #time5 = time()
                    obj = str(int(cdf.ix[i,'OBJID']))
                    srch = os.path.join(stampdir,'srch'+obj+'.gif')
                    temp = os.path.join(stampdir,'temp'+obj+'.gif')
                    diff = os.path.join(stampdir,'diff'+obj+'.gif')
                    #time6 = time()
                    
                    srcharray[i] = srch
                    #rdf.ix[i,'srchstamp'] = srch
                    #rdf.ix[i,'srchstamp'] = 'NOSTAMP'
                    
                    temparray[i] = temp
                    #rdf.ix[i,'tempstamp'] = temp
                    #rdf.ix[i,'tempstamp'] = 'NOSTAMP'
                    
                    diffarray[i] = diff
                    #rdf.ix[i,'diffstamp'] = diff
                    #rdf.ix[i,'diffstamp'] = 'NOSTAMP'
                    #time7 = time()
                    #print '1-2',time2-time1
                    #print '2-3',time3-time2
                    #print '3-4',time4-time3
                    #print '4-5',time5-time4
                    #print '5-6',time6-time5
                    #print '6-7',time7-time6
                    #sys.exit()
    
    rdf['srchstamp'] = srcharray
    rdf['tempstamp'] = temparray
    rdf['diffstamp'] = diffarray

    rdf['cutflag'] = rdf['cutflag'].astype(int)

    spl = fitsname.split('.fits')
    newfits = spl[0]+'stamps.fits'

    newfile = fitsio.FITS(newfits,'rw')
    newfile.write(rdf.to_records(index=False),clobber=True)
    newfile.close()
                    
#    for i in list(rdf.index.values):
#        nite = str(rdf.ix[i,'NITE'])
#        exp = str(rdf.ix[i,'EXPNUM'])
#        dp = 'dp'+str(season)
#        band = str(rdf.ix[i,'BAND'])
#        ccd = '%02d' % rdf.ix[i,'CCDNUM']
#        stampstar = os.path.join(expdir,nite+'/'+exp+'/'+dp+'/'+band+'_'+ccd+'/'+'stamps_*')
        

# f1= open(str(outdir)+'/'+'allcandidates.txt', 'w')
# header1 = 'SNID, ' + ' RA, ' + ' DEC, ' + ' CandType,' +  ' NumEpochs, ' + ' NumEpochsml, ' + ' LatestNiteml' 
# f1.write(header1)
# for i in range(0,numofcan):
#     Cand =(reals.data.SNID == urID[i])
#     savedata(reals,urID[i],outdir,triggerid)
# #    if not debug:
# #        continue
#     line = str(urID[i]) + ", " + str(reals.data.RA[Cand][1]) + ", " + str(reals.data.DEC[Cand][1]) + ", " + str(reals.data.CandType[Cand][1]) + ", " + str(reals.data.NumEpochs[Cand][1]) + ", " + str(reals.data.NumEpochsml[Cand][1]) + ", " + str(reals.data.LatestNiteml[Cand][1]) + "\n"
#     table1 = np.array([[int(urID[i]),reals.data.RA[Cand][1],reals.data.DEC[Cand][1],int(reals.data.CandType[Cand][1]),int(reals.data.NumEpochs[Cand][1]),int(reals.data.NumEpochsml[Cand][1]),int(reals.data.LatestNiteml[Cand][1])]])
#     print table1
#     f1.write(line)
#     filename = 'Candidate_'+str(int(urID[i])) + '.txt'
#     htmlfilename = 'Candidate_'+str(int(urID[i])) + '.html'
#     header ='BAND ' + 'x ' + 'y ' + 'Mag ' + 'Nite ' + 'MJD ' + 'Season ' + 'Object ' + 'Exposure ' + 'Field ' + 'CCDNUM'
#     nobs = len(reals.data.BAND[Cand])
#     seasoncol = np.ones((nobs,), dtype = np.int)*int(season)
#     print seasoncol
#     table = np.column_stack((reals.data.BAND[Cand], reals.data.XPIX[Cand], reals.data.YPIX[Cand], reals.data.MAG[Cand], reals.data.OBSNITE[Cand], reals.data.MJD[Cand], seasoncol))
#     np.savetxt(str(outdir)+'/'+filename, table, fmt = '%s', header = header)
#     htmlcode = HTML.table(table.tolist(),header_row = header.split(' '))
#     htmlcode1 = HTML.table(table1.tolist(), header_row = header1.split(', '))
#     f = open(str(outdir)+ '/'+  htmlfilename, 'w')
#     f.write(htmlcode1)
#     f.write(htmlcode)
#     #Collect Stamps for observations of this candidate#
#     thiscand_stampsdir = outstamps  + '/' + str(int(urID[i]))
#     if not os.path.exists(thiscand_stampsdir):
#         os.mkdir(thiscand_stampsdir)
# #Create Stamps_table with 8 columns and nobs rows all empty#
#     stampstable = ([[None]*8])
#     stampsheader = 'Filter, ' + 'Object ID, ' + 'Nite, ' + 'MJD, ' + 'Search, ' + 'Template, ' + 'Difference, ' + 'AutoScan Score,'
#     for j in range(0,nobs):
#         thisobs_nite = str(int(realss.OBSNITE[Cand][j]))
#         thisobs_band = realss.BAND[Cand][j]
#         ccdnum = int(realss.OBSCCDNUM[Cand][j])
#         objexpnum = int(realss.EXPNUM[Cand][j])
#         if int(realss.OBSCCDNUM[Cand][j]) <10:
#             ccdnum = '0' + str(int(realss.OBSCCDNUM[Cand][j]))
# #        expdir = "/data/des41.a/data/marcelle/diffimg/local-runs"
#         thisobs_ID = realss.OBJID[Cand][j]
#         a = expdir + '/' + thisobs_nite + '/'+ str(objexpnum)+'/dp' + str(season) + '/' + thisobs_band + '_' + str(ccdnum) + '/stamps*'
#         print a
#         thisobs_stampsdir = glob.glob(a)[0]
#         print thisobs_stampsdir
# ### MUST UNDO THESE COMMENTS ONCE OBSID IS DEFINED TERM ###
#         filenamediff = thisobs_stampsdir + '/diff' + str(thisobs_ID) + '.gif'
#         filenamesrch = thisobs_stampsdir + '/srch' + str(thisobs_ID) + '.gif'
#         filenametemp = thisobs_stampsdir + '/temp' + str(thisobs_ID) + '.gif' 
#         if thisobs_ID != 0:
#             shutil.copy(filenamediff, thiscand_stampsdir)
#             shutil.copy(filenamesrch, thiscand_stampsdir)
#             shutil.copy(filenametemp, thiscand_stampsdir)
#             path1 = thiscand_stampsdir + '/srch' + str(thisobs_ID) + '.gif'
#             print path1
#             search= image('', 'stamps/' + str(int(urID[i]))  + '/srch' + str(thisobs_ID) + '.gif')
#             temp  = image('', 'stamps/' + str(int(urID[i])) + '/temp' + str(thisobs_ID) + '.gif')
#             diff  = image('', 'stamps/' + str(int(urID[i])) + '/diff' + str(thisobs_ID) + '.gif')
            
#         if thisobs_ID == 0:
#             search = 'no search'
#             temp = 'no temp'
#             diff = 'no diff'
# #Replace in the empty spaces in table with values/pictures#
#         stampstable[j][0] = realss.BAND[Cand][j]
#         stampstable[j][1] = realss.OBJID[Cand][j]
#         stampstable[j][2] = realss.OBSNITE[Cand][j]
#         stampstable[j][3] = realss.MJD[Cand][j]
#         stampstable[j][4] = search
#         stampstable[j][5] = temp
#         stampstable[j][6] = diff
#         stampstable[j][7] = realss.PHOTPROB[Cand][j]
#         stampstable.append([None] * 8)
        
#     htmlcode2 = HTML.table(stampstable, header_row= stampsheader.split(', '))
#     f.write(htmlcode2)
#     #Making Plot of Flux vs MJD for each Candidate#
#     fvmplottable = ([[None]*len(ubands)])
#     for b in range(0, len(ubands)):
#         Bandcand = realss.BAND[Cand]
#         Band = Bandcand == ubands[b]
#         Flux = realss.FLUXCAL[Cand][Band]
#         MJD = realss.MJD[Cand][Band]
#         Fluxerr = realss.FLUXCALERR[Cand][Band]
#         plt.scatter(MJD,Flux, color = colors[b])
#         plt.errorbar(MJD,Flux, yerr=Fluxerr,color=colors[b], ls = 'none')
#         plt.xlabel('MJD')
#         plt.ylabel('Flux')
#     plt.title('Flux vs. MJD for candidate'  + str(int(urID[i]))+ ' iband=r, zband=g')
#     plt.savefig(outdir +'/plots/lightcurves/FluxvsMJD_for_cand_' + str(int(urID[i]))+'.png')
#     plt.clf()
#     plotpic = image('', 'plots/lightcurves/FluxvsMJD_for_cand_' + str(int(urID[i])) + '.png')
#     fvmplottable[0][b] = plotpic
#     htmlcode3 = HTML.table(fvmplottable)
#     f.write(htmlcode3)

                     
# f1.close()  
    
# print "SUCCESS"    


             


# #Save output files in specific and organized directories#

# ###Setup a display of the search, template and difference images for each candidate. Automatically save for top n candidates. Give option to user to pick candidate by ID to display###
