import glob
import os
import pyncl.nclfuncs as ncl

def run():
    print('starting: namelist.py')   # line 13 NCL script
    o = os.getenv('OBS')
    case_sens = os.getenv('MACHINE')

    if o:
        obsflag = True
    else:
        obsflag = False

    if not os.path.exists('namelist_byvar/'):
        os.mkdir('namelist_byvar')

    # original NCL script goes to some trouble to deal with blank lines and to
    # determine number of simulations so that it can allocate some arrays of that size.
    # we're going to simplify that a bit here.
    na = []
    f = open('namelist')       # line 28 NCL script
    for line in f:
        if line != "":
            na.append(line)
    f.close()
    nsims = len(na)

    os.putenv('NSIM', str(nsims))      # line 40 NCL script

    names = []
    paths = []
    syear = []
    eyear = []

    delim = '|'                # line 48 NCL script
    for gg in na:
        tokens = gg.split(delim)
        names.append(tokens[0].strip())
        paths.append(tokens[1].strip())
        syear.append(tokens[2].strip())
        eyear.append(tokens[3].strip())

    # if path ends in.nc remove it. (It will get appended to the end of the path
    # automatically when searching below.)
    for gg in range(len(paths)):              # line 59 NCL script
        paths[gg] = paths[gg].replace('.nc','')

    #----- Read in namelist_obs, and check number of supplied Observational datasets ------

    maxnumobs = 0    # set maximum number of obs datasets per variable. if(obsflag).eq.True, this will likely get altered.

    if obsflag:      # line 69 NCL script
        # again, we're going to simplify some of the reading of the namelist_obs file...
        na = []
        f = open('namelist_obs')       # line 78 NCL script
        for line in f:
            line = line.strip()
            if len(line) > 0:
                na.append(line)
        f.close()
        nobs = len(na)

        vnamesB = []
        namesB = []
        pathsB = []
        syearB = []
        eyearB = []

        strFillVal = ncl.default_fillvalues['string']
        intFillVal = ncl.default_fillvalues['int32']
        for gg in na:                # line 83 in NCL script
            tokens = gg.split(delim)
            vnamesB.append( tokens[0].strip())
            t = tokens[1].strip()
            namesB.append( (strFillVal, t)[t != ''] )
            t = tokens[2].strip()
            pathsB.append( (strFillVal, t)[t != ''] )
            t = tokens[3].strip()
            syearB.append( (strFillVal, t)[t != '']  )
            t = tokens[4].strip()
            eyearB.append( (strFillVal, t)[t != '']  )

        maxnumobs = max(vnamesB.count('TS'), vnamesB.count('PSL'), vnamesB.count('TREFHT'),
                        vnamesB.count('PRECT'), vnamesB.count('MOC'), vnamesB.count('SNOWDP'),
                        vnamesB.count('aice_nh'), vnamesB.count('aice_sh'))

        #  check to see if any names are duplicated. If they are, add a "_2", "_3" to the name
        # this is needed so that each output .nc file has a different name
        for gg in namesB:           # line 102 NCL script
            if namesB.count(gg) > 1:
                count = 0
                for i in range(len(namesB)):
                    if namesB[i] == gg:
                        count += 1
                        if count > 1:
                            namesB[i] = namesB[i] + '_' + str(count-1)

        f = open('obs_maxnum', 'w')       # line 113 NCL script
        f.write(str(maxnumobs))
        f.close()

    # ----- TS section---------------
    namelist_ts = []
    if obsflag:
        ts_i = indEQ(vnamesB, "TS")
        if len(ts_i) > 0:
            for gg in ts_i:
                namelist_ts.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            # fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
            suffix = 1
            while len(namelist_ts) < maxnumobs:
                namelist_ts.append(f'{namesB[ts_i[0]]}_{str(suffix)} | {pathsB[ts_i[0]]} | {syearB[ts_i[0]]} | {eyearB[ts_i[0]]}')
                suffix += 1

            f = open('obs_ts', 'w')         # line 135 NCL script
            f.write(namelist_ts[0] + '\n')
            f.close()
        else:
            namelist_ts = ["missing"]*maxnumobs

    if case_sens:
        tstring0 = "TS_,.TS.,ts_,.ts.,t_surf_,t_surf.,sst.,sst_"
        tstring1 = "TS,ts,t_surf,sst"
    else:
        tstring0 = "TS_,.TS.,t_surf_,t_surf.,sst.,sst_"
        tstring1 = "TS,t_surf,sst"

    for gg,path in enumerate(paths):         # line 149 NCL script        
        fsst = myGlob(paths[gg]+'*{'+tstring0+'}*.nc')
        if len(fsst) == 1:
            namelist_ts.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = path.replace("/*/", "/{"+tstring1+"}/")     # explicitly specify TS,ts in directory structure to eliminate "/tsmin/" being used
            namelist_ts.append( f'{names[gg]} | {tpath}{{{tstring0}}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_ts', 'w')   # line 161 NCL script
    for name in namelist_ts:
        f.write(name + '\n')
    f.close()

    #------- PSL section----------------------------  
    namelist_psl = []          # line 163 NCL script
    if obsflag:   
        psl_i = indEQ(vnamesB, "PSL")
        if len(psl_i) > 0:
            for gg in psl_i:
                namelist_psl.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            # fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
            suffix = 1
            while len(namelist_psl) < maxnumobs:
                namelist_psl.append(f'{namesB[psl_i[0]]}_{str(suffix)} | {pathsB[psl_i[0]]} | {syearB[psl_i[0]]} | {eyearB[psl_i[0]]}')
                suffix += 1

            f = open('obs_psl', 'w')
            f.write(namelist_psl[0] + '\n')
            f.close()
        else:
            namelist_psl = ['missing']*maxnumobs

    if case_sens:           # line 186 NCL script
        tstring0 = "PSL_,PSL.,psl_,psl.,slp.,slp_"
        tstring1 = "PSL,psl,SLP,slp"
    else:
        tstring0 = "PSL_,PSL.,slp.,slp_"
        tstring1 = "PSL,slp"

    for gg, path in enumerate(paths):
        fsst = myGlob(path+'*{'+tstring0+'}*.nc')
        if len(fsst) == 1:
            namelist_psl.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{"+tstring1+"}/")
            namelist_psl.append(f'{names[gg]} | {tpath}*{{{tstring0}}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_psl', 'w')       # line 205 NCL script
    for line in namelist_psl:
        f.write(line + '\n')
    f.close()

    # ------- TREFHT section----------------------------
    namelist_trefht = []            # line 207 NCL script
    if obsflag:
        trefht_i = indEQ(vnamesB, 'TREFHT')
        if len(trefht_i) > 0:
            for gg in trefht_i:
                namelist_trefht.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            # fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
            suffix = 1
            while len(namelist_trefht) < maxnumobs:
                namelist_trefht.append(f'{namesB[trefht_i[0]]}_{str(suffix)} | {pathsB[trefht_i[0]]} | {syearB[trefht_i[0]]} | {eyearB[trefht_i[0]]}')
                suffix += 1

            f = open('obs_trefht', 'w')     # line 224 NCL script
            f.write(namelist_trefht[0] + '\n')
            f.close()
        else:
            namelist_trefht = ['missing']*maxnumobs

    if case_sens:         # line 230 NCL script
        tstring0 = "TREFHT_,TREFHT.,tas.,tas_,t_ref.,t_ref_,T2.,T2_"
        tstring1 = "TREFHT,tas,t_ref,T2"
    else:
        tstring0 = "TREFHT_,TREFHT.,tas.,tas_,t_ref.,t_ref_,T2.,T2_"
        tstring1 = "TREFHT,tas,t_ref,T2"

    for gg,path in enumerate(paths):
        fsst = myGlob(paths[gg]+"*{"+tstring0+"}*.nc")
        if len(fsst) == 1:
            namelist_trefht.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{"+tstring1+"}/")
            namelist_trefht.append(f'{names[gg]} | {tpath}*{{{tstring0}}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_trefht', 'w')    # line 224 NCL script
    for line in namelist_trefht:
        f.write(line + '\n')
    f.close()

    # ------- PRECT section--(more complicated due to PRECC+PRECL, + pr being a common 2 letter combination)------
    namelist_prect = []      # line 250 NCL script
    if obsflag:
        prect_i = indEQ(vnamesB, 'PRECT')
        if len(prect_i) > 0:
            for gg in prect_i:
                namelist_prect.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            # fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-ouprect
            suffix = 1
            while len(namelist_prect) < maxnumobs:
                namelist_prect.append(f'{namesB[prect_i[0]]}_{str(suffix)} | {pathsB[prect_i[0]]} | {syearB[prect_i[0]]} | {eyearB[prect_i[0]]}')
                suffix += 1

            f = open('obs_prect', 'w')  # line 224 NCL script
            f.write(namelist_prect[0] + '\n')
            f.close()
        else:
            namelist_prect = ['missing']*maxnumobs

    pstring = ["pr_*", "pr.*", "_pr_*", ".pr.*", "PRECT.*", "PRECT_*", "PRECC.*", "PRECC_*", "precip_*", "precip.*",
               "prcp_*", "prcp.*", "prate_*", "prate.*"]
    for gg in range(nsims):     # line 272 NCL script
        #RLB the NCL script employs a different means for missing paths. It relies on
        #RLB that the systemfunc"ls ...") command will return a 1 element string result
        #RLB set to "missing" if the path is not matched. For every variable in pstring
        #RLB that is not matched, it (re)writes "missing" into the "current" position
        #RLB in the namelist_prect array until a match (if any) is found.
        #RLB Here, we are using a list and appending results -- there is no "current"
        #RLB index into namelist_prect. So we'll check at the end of the inner loop if
        #RLB anything has been found for the simulation, and force a "missing" into the
        #RLB list if not.
        foundPrect = False
        for pstringhh in pstring:
            fsst = myGlob(paths[gg]+"*"+pstringhh+".nc")
            if len(fsst) == 1:
                if pstringhh == 'PRECC.*' or pstringhh == 'PRECC_*':
                    tpath = paths[gg].replace("/*/", "/{PRECC,PRECL}/")
                    namelist_prect.append(f'{names[gg]} | {tpath}*{{PRECC,PRECL}}*.nc | {syear[gg]} | {eyear[gg]}')
                else:
                    namelist_prect.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
                foundPrect = True
                break
            else:
                if pstringhh == 'PRECC.*' or pstringhh == 'PRECC_*':       # line 293 NCL script
                    tpath = paths[gg].replace("/*/", "/{PRECC,PRECL}/")
                    namelist_prect.append(f'{names[gg]} | {tpath}*{{PRECC,PRECL}}*.nc | {syear[gg]} | {eyear[gg]}')
                else:
                    tpath = ''
                    if pstringhh == 'pr_*' or pstringhh == 'pr.*' or pstringhh == '*_pr_*' or pstringhh == '*.pr.*':
                        tpath = paths[gg].replace("/*/", "/pr/")
                    if pstringhh == 'PRECT.*' or pstringhh == 'PRECT_*':
                        tpath = paths[gg].replace("/*/", "/PRECC/")
                    if len(tpath) > 0:
                        namelist_prect.append(f'[{names[gg]} | {tpath}*{pstringhh}*.nc | {syear[gg]} | {eyear[gg]}')
                    else:
                        namelist_prect.append(f'{names[gg]} | missing')
                foundPrect = True
                break
        if not foundPrect:
            namelist_prect.append(f'{names[gg]} | missing | {syear[gg]} | {eyear[gg]}')   # file is missing..

    f = open('namelist_byvar/namelist_prect', 'w')     # line 318 NCL script
    for line in namelist_prect:
        f.write(line + '\n')
    f.close()

    # ----- SNOWDP section---------------
    namelist_snowdp = []     # line 356 NCL script
    if obsflag:
        snowdp_i = indEQ(vnamesB, "SNOWDP")
        if len(snowdp_i) > 0:
            for gg in snowdp_i:
                namelist_snowdp.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            suffix = 1
            while len(namelist_snowdp) < maxnumobs:
                namelist_snowdp.append(f'{namesB[snowdp_i[0]]}_{str(suffix)} | {pathsB[snowdp_i[0]]} | {syearB[snowdp_i[0]]} | {eyearB[snowdp_i[0]]}')
                suffix += 1

            f = open('obs_snowdp', 'w')       # line 337 NCL script
            f.write(namelist_snowdp[0] + '\n')
            f.close()
        else:
            namelist_snowdp = ['missing']*maxnumobs

    sn_string = 'SNOWDP_,SNOWDP.,snd_,snd.'
    for gg in range(nsims):
        fsst = myGlob(paths[gg]+'*{'+sn_string+'}*.nc')
        if len(fsst) == 1:
            namelist_snowdp.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{SNOWDP,snd}/")     # explicitly specify SNOWDP/snd in directory structure to eliminate "/sndmin/" being used
            namelist_snowdp.append(f'{names[gg]} | {tpath}*{{SNOWDP_,SNOWDP.,snd_,snd.}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_snowdp', 'w')    # line 354 NCL script
    for line in namelist_snowdp:
        f.write(line + '\n')
    f.close()

    # ------- MOC section----------------------------
    namelist_moc = []        # line 356 NCL script
    if obsflag:
        moc_i = indEQ(vnamesB, 'MOC')
        if len(moc_i):
            for gg in moc_i:
                namelist_moc.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            suffix = 1
            while len(namelist_moc) < maxnumobs:
                namelist_moc.append(f'{namesB[moc_i[0]]}_{str(suffix)} | {pathsB[moc_i[0]]} | {syearB[moc_i[0]]} | {eyearB[moc_i[0]]}')
                suffix += 1

            f = open('obs_moc', 'w')
            f.write(namelist_moc[0] + '\n')
            f.close()
        else:
            namelist_moc = ['missing']*maxnumobs

    for gg in range(nsims):
        fsst = myGlob(paths[gg]+"*{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}*.nc")   # /dev/null suppresses all standard error output
        if len(fsst) == 1:
            namelist_moc.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{MOC,msftmyz,stfmmc}/")
            namelist_moc.append(f'{names[gg]} | {tpath}*{{MOC_,MOC.,msftmyz.,msftmyz_,stfmmc.,stfmmc_}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_moc', 'w')      # line 390 NCL script
    for line in namelist_moc:
        f.write(line + '\n')
    f.close()

    # ------- aice_nh section----------------------------
    namelist_aice_nh = []       # line 393 NCL script
    if obsflag:
        aice_nh_i = indEQ(vnamesB, "aice_nh") + indEQ(vnamesB, "AICE_NH")
        if len(aice_nh_i) > 0:
            for gg in aice_nh_i:
                namelist_aice_nh.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            suffix = 1
            while len(namelist_aice_nh) < maxnumobs:
                namelist_aice_nh.append(f'{namesB[aice_nh_i[0]]}_{str(suffix)} | {pathsB[aice_nh_i[0]]} | {syearB[aice_nh_i[0]]} | {eyearB[aice_nh_i[0]]}')
                suffix +=1

            f = open('obs_aice_nh', 'w')
            f.write(namelist_aice_nh[0] + '\n')
            f.close()
        else:
            namelist_aice_nh = ['missing']*maxnumobs

    for gg in range(nsims):
        fsst = myGlob(paths[gg]+'*{aice_nh.,aice.,sic_,sic.,.CN.,_CN_}*.nc')
        if len(fsst) == 1:
            namelist_aice_nh.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{aice,sic,aice_nh,CN}/")
            namelist_aice_nh.append(f'{names[gg]} | {tpath}*{{aice_nh.,aice.,sic_,sic.,.CN.,_CN_}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_aice_nh', 'w')    # line 427 NCL script
    for line in namelist_aice_nh:
        f.write(line + '\n')
    f.close()

    # ------- aice_sh section----------------------------
    namelist_aice_sh = []         # line 429 NCL script
    if obsflag:
        aice_sh_i = indEQ(vnamesB, 'aice_sh') + indEQ(vnamesB, 'AICE_SH')
        if len(aice_sh_i) > 0:
            for gg in aice_sh_i:
                namelist_aice_sh.append(f'{namesB[gg]} | {pathsB[gg]} | {syearB[gg]} | {eyearB[gg]}')
            suffix = 1
            while len(namelist_aice_sh) < maxnumobs:     #fill in the missing obs rows with the first obs file, altering the name slightly for .nc write-outs
                namelist_aice_sh.append(f'{namesB[aice_sh_i[0]]}_{str(suffix)} | {pathsB[aice_sh_i[0]]} | {syearB[aice_sh_i[0]]} | {eyearB[aice_sh_i[0]]}')
                suffix += 1

            f = open('obs_aice_sh', 'w')   # line 446 NCL script
            f.write(namelist_aice_sh[0] + '\n')
            f.close()
        else:
            namelist_aice_sh = ['missing']*maxnumobs

    for gg in range(nsims):
        fsst = myGlob(paths[gg]+"*{aice_sh.,aice.,sic_,sic.,.CN.,_CN_}*.nc")
        if len(fsst) == 1:
            namelist_aice_sh.append(f'{names[gg]} | {fsst[0]} | {syear[gg]} | {eyear[gg]}')   # grab first file
        else:
            tpath = paths[gg].replace("/*/", "/{aice,sic,aice_sh,CN}/")
            namelist_aice_sh.append(f'{names[gg]} | {tpath}*{{aice_sh.,aice.,sic_,sic.,.CN.,_CN_}}*.nc | {syear[gg]} | {eyear[gg]}')

    f = open('namelist_byvar/namelist_aice_sh', 'w')   # line 463 NCL script
    for line in namelist_aice_sh:
        f.write(line + '\n')
    f.close()

    # ----------------------------------------------------------------------------
    print('Finished: namelist.py')


# -----------------------------------------------------------------------------------
# Emulates NCL's ind( .eq. ) function, when handed a set object
#
def indEQ(collection, item):
    indices = []
    for index, c in enumerate(collection):
        if c == item:
            indices.append(index)
    return indices

# -----------------------------------------------------------------------------------
# Implementation of "myGlob()"
#
# Python's glob.glob() function does not support so-called "brace-expansion", which is
# utilized heavily by this script. myGlob() provides the required behavior. It makes use of
# two functions, 'getitem()' and 'getgroup(), which have been taken from the Python example
# at:
#
#   https://rosettacode.org/wiki/Brace_expansion#Python
#
# Note that the rosettacode.org website states that: "Content is available under GNU Free Documentation License 1.2
# unless otherwise noted.":
#
#   https://rosettacode.org/wiki/Rosetta_Code (bottom of page, as of 9/20/2018)
#
# Beyond rosettacode.org, there is no author to directly cite.
#
def getitem(s, depth=0):
    out = [""]
    while s:
        c = s[0]
        if depth and (c == ',' or c == '}'):
            return out, s
        if c == '{':
            x = getgroup(s[1:], depth + 1)
            if x:
                out, s = [a + b for a in out for b in x[0]], x[1]
                continue
        if c == '\\' and len(s) > 1:
            s, c = s[1:], c + s[1]

        out, s = [a + c for a in out], s[1:]

    return out, s


def getgroup(s, depth):
    out, comma = [], False
    while s:
        g, s = getitem(s, depth)
        if not s: break
        out += g

        if s[0] == '}':
            if comma: return out, s[1:]
            return ['{' + a + '}' for a in out], s[1:]

        if s[0] == ',':
            comma, s = True, s[1:]

    return None

def myGlob(path):
    for s in path.split('\n'):
        expandedPath = getitem(s)[0]

    filelist = []
    for file in expandedPath:
        filelist += glob.glob(file)

    if len(filelist) == 0:  # this emulates the systemfunc returned value when no match is found
        filelist.append('missing')

    return filelist

# -----------------------------------------------------------------------------------

if __name__ == '__main__':
    run()