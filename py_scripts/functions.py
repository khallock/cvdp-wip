import pyncl.nclfuncs as ncl
import numpy as np
import xarray as xr
import glob
import math
import os

def create_empty_array( yS, yE, mS, mE, opt_type):
    '''create blank array for use when something may be/is wrong.'''
    if yS is None or yE is None:
        yS.values[:] = 1
        yE.values[:] = 50

    timeT = ncl.yyyymm_time(yS, yE, int)
    time = timeT.sel(time=slice(yS*100+mS, yE*100+mE))

    if opt_type == 'time_lat_lon':
        blank_array_values = np.empty((time.sizes['time'], 90, 180), dtype=np.float32)
        lat = xr.DataArray(np.linspace(-89,89,90), dims=('lat'), attrs={'units':'degrees_north'})
        lon = xr.DataArray(np.linspace(0,358,180), dims=('lon'), attrs={'units':'degrees_east'})
        blank_array = xr.DataArray(blank_array_values, dims=('time', 'lat', 'lon'),
                                   coords={'time':time, 'lat':lat, 'lon':lon})

    elif opt_type == 'time_lev_lat':
        blank_array_values = np.empty((time.sizes['time'], 41, 90), dtype=np.float32)
        lat = xr.DataArray(np.linspace(-89,89,90), dims=('lat'), attrs={'units':'degrees_north'})
        lev = xr.DataArray(np.linspace(0,500,41), dims=('lev'), attrs={'units':'m', 'positive':'down'})
        blank_array = xr.DataArray(blank_array_values, dims=('time', 'lev', 'lat'),
                                   coords={'time' : time, 'lev' : lev, 'lat' : lat})

    blank_array.attrs['units'] = ''
    blank_array.attrs['is_all_missing'] = True

    return blank_array

def data_read_in(zpath, vn, yearS, yearE):
    '''
    read in atmospheric / land data from selected files
    assign time coordinate variables, check for issues with the array, assign _FillValue ( if needed)
    assign dimension names(for ease - of - use), check and modify units

    vname settings at top of this script can be modified if a different variable name is
    encountered. For instance, if a TS data file has the TS array named as "sfc_t", one
    could add "sfc_t" to the vname TS coding as follows:
        if (vn.eq."TS") then
            vname = (/ "TS", "ts", "sst", "sfc_t" /)
        end if
    '''
    # path for TS file(s), variable name, start year, and end year are read in.

    arr = None
    vname = None

    if vn == 'TS':
        vname = ("TS","ts","sst","t_surf","skt")
    elif vn == "PSL":
        vname = ("PSL","psl","slp","SLP","prmsl","msl","slp_dyn")
    elif vn == "TREFHT":
        vname = ("TREFHT","tas","temp","air","temperature_anomaly","temperature","t2m","t_ref","T2","tempanomaly")
    elif vn == "PRECT":
        vname = ("PRECC","PRECL","PRECT","pr","PPT","ppt","p","P","precip","PRECIP","tp","prcp","prate")
    elif vn == "SNOWDP":
        vname = ("SNOWDP","snd")

    if zpath is None:
        print(f'File missing, creating blank array of data. View {vn} namelist for details.')
        arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        if '*' in zpath or '{' in zpath:     # check for "*" and "{" denoting multiple files
            tfiles = glob.glob(zpath)
            if vn == 'PRECT':  # special section for precip, as might need to do PRECC+PRECL
                b = xr.open_dataset(tfiles[0])
                if 'PRECC' in b or 'PRECL' in b:    # PRECC/PRECL section
                    fils_precc = []
                    fils_precl = []
                    for f in b:
                        if 'PRECC' in b:
                            fils_precc.append(f)
                        elif 'PRECL' in b:
                            fils_precl.append(b)

                    if len(fils_precc) == 0 or len(fils_precl) == 0:
                        print('Fatal: Need both PRECC and PRECL file(s), creating blank array')
                        print(fils_precc)
                        print(fils_precl)
                        arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
                        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
                        eydata = yearE     # as data array is totally missing..
                        smdata = 1
                        emdata = 12
                        ##
                        ## In the NCL code, there is a break here, however I don't
                        ## think it works as intended. Replacing with a return  -- RLB
                        ##
                        return

                    with xr.open_mfdataset(fils_precc) as c:
                        arr1 = c.data_vars['PRECC']
                    with xr.open_mfdataset('fils_precl') as c:
                        arr2 = c.data_vars['PRECL']
                    arr = xr.concat(arr1, arr2)
                    arr.attrs['long_name'] = 'Large-scale (stable) + convective precipitation rate (liq + ice)'

                else:   # pr, ppt, PPT, PRECT multiple/single file read-in here..
                    c = xr.open_mfdataset(tfiles)
                    for v in vname:
                        if v in b:
                            arr = c.data_vars[v]
                            break
            else:
                c = xr.open_mfdatatset(tfiles)
                for v in vname:
                    if v in c:
                        arr = c.data_vars[v]
                        break

            cpathS = tfiles[0]
            cpathE = tfiles[len(tfiles)-1]
            sydata = cpathS[len(cpathS)-17:len(cpathS)-13]
            smdata = cpathS[len(cpathS)-13:len(cpathS)-11]
            eydata = cpathE[len(cpathE)-10:len(cpathE)-6]
            emdata = cpathE[len(cpathE)-6:len(cpathS)-4]

        else:    # single file case...
            c = xr.open_dataset(zpath)
            for v in vname:
                if v in c:
                    arr = c.data_vars[v]
                    break
            sydata = zpath[len(zpath)-17:len(zpath)-13]
            smdata = zpath[len(zpath)-13:len(zpath)-11]
            eydata = zpath[len(zpath)-10:len(zpath)-6]
            emdata = zpath[len(zpath)-6:len(zpath)-4]

    if arr is None:
        print(f'Variable {vn} not found. Examine input file {zpath}. Creating empty array and continuing')
        arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    if arr.values.dtype == 'int16':
        arr.values = arr.values.astype('float32')

    #### HOW TO DEAL WITH MISSING_VALUES IN XARRAY?  --RLB
    if '_FillValue' not in arr.attrs:   # assign _FillValue if one is not present
        if 'missing_value' in arr.attrs:
            arr.attrs['_FillValue'] = arr.attrs['missing_value']
        else:
            arr.attrs['_FillValue'] = ncl.default_fillvalues[arr.dtype.name]

    if [True for dimsize in arr.sizes if dimsize == 1]:
        arr = arr.squeeze()

    if len(arr.dims) <= 2:
        print('Possible curvilinear (or unstructured) grid detected. The CVDP cannot analyze curvilinear data. Please regrid to a rectilinear grid for inclusion in CVDP comparisons.')
        print(f'Input file: {zpath}')
        print('Setting array to all missing')
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12

    ###RLB - is this the correct idiom in Xarrays for NCL's dimension (re)naming?
    arr.rename({arr.dims[0] : 'time'})
    arr.rename({arr.dims[1] : 'lat'})
    arr.rename({arr.dims[2] : 'lon'})

    if 'valid_range' in arr.attrs:  # check to make sure data is in valid range. Reset to stay within the valid range if needed.
        arr.data = xr.where(arr.data < arr.attrs['valid_range'][0], arr.attrs['valid_range'][0], arr.data)
        arr.data = xr.where(arr.data > arr.attrs['valid_range'][1], arr.attrs['valid_range'][1], arr.data)

    ###RLB - what are we doing with _FillValues?
    ###  if arr.abs() > 1.e20:   # check for inf values or values way out of range, reset to _FillValue.
    ###      print(f'Values greater than 1.e20 or less than -1.e20 detected in {zpath}, resetting to _FillValue')
    ###      arr = xr.where(arr.abs() > 1.e20, arr@_FillValue,arr)
    print(arr)


    if yearS < int(sydata) or yearE > int(eydata):
        print(f'Requested {yearS}-{yearE} time span is outside the input file {zpath} time span of {sydata}-{eydata}')
        print('Setting array to all missing')
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        timeT =  ncl.yyyymm_time(sydata, eydata, "integer")
        time = timeT[slice(sydata*100+smdata, eydata*100+emdata)]
        if 'time' in arr.coords:
            arr.drop('time')
        dimz = arr.values.shape
        if dimz[0] == time.values.shape[0]:
            #### THIS MAY NOT BE RIGHT -- MAY NEED TO REASSIGN ARR RESULTS OF ASSIGN_COORDS
            arr.assign_coords(time=time)
        else:
            print('Possible mismatch detected between time specified in file name and file variables, setting array to missing')
            print(f'File = {zpath}')
            print(f'Read from file name: {min(time)}-{max(time)}')
            arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
            sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
            eydata = yearE     # as data array is totally missing..
            smdata = 1
            emdata = 12

    if arr.coords['lat'][0] >= 0:
        farr = arr[slice(yearS*100+1, yearE*100+12), ::-1, :]   # flip the latitudes
    else:
        farr = arr[slice(yearS*100+1, yearE*100+12), :, :]

    mocheck = np.array( (yearS*100+1)-min(farr.coords['time'], (yearE*100+12) - max(farr.coords[time])) )
    if [True for mon in mocheck if mon != 0]:    # previously: if (mod(dimsizes(farr&time),12).ne.0) then
        if mocheck[0] != 0:
            print("First requested year is incomplete")
        if mocheck[1] != 0:
            print("Last requested year is incomplete")
        print(f'Incomplete data year(s) requested for file {zpath}, printing out time and creating blank array')
        print(f'Time requested: {yearS}-{yearE}')
        print(farr.coords['time'])
        farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    if farr.coords['lon'][0] < 0:
        farr = ncl.lonFlip(farr)     # lon flip    #####RLB: THIS IS STUBBED FOR RIGHT NOW!!!

    if farr.coords['lon'].min() < 0 or farr.coords['lon'].max() > 360:
        print(farr.coords['lon'])
        print(f'path = {zpath}')
        print("Fatal: Longitudes not in expected 0-360E range, creating blank array")
        farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    if vn == 'TREFHT' or vn == 'TS':      # units check
        if farr.attrs['units']  == 'K' or farr.attrs['units'] == 'Kelvin' or \
                farr.attrs['units'] == 'deg_k' or farr.attrs['units'] == 'deg_K':
            if farr.max() >= 100:    # data sets can be anomalies with units of K, so check for range before subtracting
                farr.values = farr.values - 273.15
            farr.attrs['units'] = 'C'
        if farr.attrs['units'] == 'degrees_C' or farr.attrs['units'] == 'degrees C' or \
                farr.attrs['units'] == 'degree_C' or farr.attrs['units'] == 'degree C':
            farr.attrs['units'] = 'C'

    if vn == 'PSL':
        if farr.attrs['units'] == 'Pa' or farr.attrs['units'] == 'Pascals' or farr.attrs['units'] == 'Pascal':
            farr.values = farr.values / 100.
            farr.attrs['units']= 'hPa'

    if vn == 'PRECT':    # convert (if necessary) to mm/day
        if farr.attrs['units'] == 'm/s' or farr.attrs['units'] == 'm s-1':
            farr.values = farr.values * 86400000.
        if farr.attrs['units'] == 'kg m-2 s-1' or farr.attrs['units'] == 'kg/m2/s' or \
                farr.attrs['units'] == 'kg/m^2/s ' or farr.attrs['units'] == 'kg/(s*m2)' or \
                farr.attrs['units'] == 'mm/s':
            farr.values = farr.values * 86400.
        if farr.attrs['units'] == 'm' or farr.attrs['units'] == 'm/month' or \
                farr.attrs['units'] == 'cm'  or farr.attrs['units'] == 'cm/month' or \
                farr.attrs['units'] == 'mm' or farr.attrs['units'] == 'mm/month':
            yr = farr.coords['time'].values.astype(int) / 100
            mo = farr.coords['time'].values.astype(int) - yr*100   ##RLB THE WHOLE EXPRESSION WAS toint()!
            days = ncl.days_in_month(yr, mo)
            for gg in range(0,farr.coords['time'].sizes-1):
                farr.values[gg,:,:] = farr.values[gg,:,:] / days(gg)  ##RLB IS THIS RIGHT? LINE 322
            if farr.attrs['units'] == 'cm' or farr.attrs['units'] == 'cm/month':
                farr.values = farr.values * 10.   # convert from cm/day to mm/day
            if farr.attrs['units'] == 'm' or farr.attrs['units'] == 'm/month':
                farr.values = farr.values * 1000.   # convert from m/day to mm/day

        if farr.attrs['units'] == 'm/day' or farr.attrs['units'] == 'm day-1':
            farr.values = farr.values * 1000.
            farr.attrs['units'] = 'mm/day'

    if vn == 'SNOWDP':
        if 'is_all_missing' not in farr.attrs:
            if farr.attrs['units'] != 'm' and farr.attrs['units'] != 'meters':
                print(f'Warning: SNOWDP/snd units may not be in meters. listed units = {farr.attrs["units"]}')

    date = farr.coords['time']         # switch time to be CF-conforming
    yyyy = date/100
    mm = date-(yyyy*100)
    days = (ncl.days_in_month(yyyy,mm))/2
    hms  = days
    hms = 0   # hours, minutes, seconds all the same (=0)
    time = ncl.cd_inv_calendar(yyyy,mm,days,hms,hms,hms,"months since "+min(yyyy)+"-01-15 00:00:00",0)
    time.attrs['long_name'] = 'Time'
    time.attrs['standard_name'] = 'time'
    time.attrs['actual_range'] = np.array(min(time),max(time))
    time.rename( {time.dims[0] : 'time'}  )
    time.coords['time'] = time
    farr.coords['time'] = time
    return farr

def data_read_in_ocean_MOC(zpath, vn, yearS, yearE):
    '''
        read in MOC ocean data from given files

        assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
        assign dimension names (for ease-of-use), check and modify units

        path for MOC file(s), variable name, start year, and end year are read in.
    '''

    arr = None
    vname = None

    if vn == 'MOC':      # line 375 in ncl version
        vname = np.array('MOC', 'msftmyz','stfmmc')

    if zpath is None:
        print(f'File missing, creating blank array of data. View {vn} namelist for details.')
        arr = create_empty_array(yearS,yearE,1,12,"time_lev_lat")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        if '*' in zpath or '{' in zpath:     # check for "*" and "{" denoting multiple files
            tfiles = glob.glob(zpath)
            c = xr.open_mfdataset(tfiles)
            for v in vname:
                if v in c:
                    dimC = c.variables['MOC'].shape
                    if v == 'MOC':    # CCSM/CESM file
                        if dimC[2] >= 2:
                            arr = c[v][:,1,:,:,:].sum(axis=1)   # select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and sum over moc_comp
                        else:
                            arr = c[v][:,1,0,:,:]               # select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and the only moc_comp dimension
                    else:                                        # CMIP file
                        arr = c[v][:,0,:,:]                     # CMIP file: 0th basin/region = atlantic_ocean (CMIP3) or atlantic_arctic_ocean (CMIP5)
                    break

            cpathS = tfiles[0]
            cpathE = tfiles[len(tfiles)-1]
            sydata = cpathS[len(cpathS)-17:len(cpathS)-14]
            smdata = cpathS[len(cpathS)-13:len(cpathS)-12]
            eydata = cpathE[len(cpathE)-10:len(cpathE)-7]
            emdata = cpathE[len(cpathE)-6:len(cpathE)-5]
        else:
            c = xr.open_dataset(zpath)
            for v in vname:
                if v in c:
                    dimC = c.variables["MOC"].shape
                    if v == 'MOC':    # CCSM/CESM file
                        if dimC[2] >= 2:
                            arr = c[v][:,1,:,:,:].sum(axis=1)  # select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region and sum over moc_comp
                        else:
                            arr = c[v][:,1,0,:,:]              # select Atl+Med+Labrador+GIN sea+Arctic+Hudson Bay transport region
                    else:                                      # CMIP file
                        arr = c[v][:,0,:,:]                    # CMIP file: 0th basin/region = atlantic_ocean (CMIP3) or atlantic_arctic_ocean (CMIP5)
                break
            sydata = zpath[len(zpath)-17:len(zpath)-14]
            smdata = zpath[len(zpath)-13:len(zpath)-12]
            eydata = zpath[len(zpath)-10:len(zpath)-7]
            emdata = zpath[len(zpath)-6:len(zpath)-5]

    if arr is None:    # line 450 in ncl version
        print(f'Variable {vn} not found. Examine input file {zpath}. Creating empty array and continuing')
        arr = create_empty_array(yearS,yearE,1,12,"time_lev_lat")

    if arr.values.dtype == 'int16':
        arr.values= arr.values.astype('float32')


    if '_FillValue' not in arr.attrs:     # assign _FillValue if one is not present
        if 'missing_value' in arr.attrs:
            arr.attrs['_FillValue'] = arr.attrs['missing_value']
        else:
            arr.attrs['_FillValue'] = ncl.default_fillvalues(arr.dtype.name)

    arr.rename( {arr.dims[0] : 'time'} )
    arr.rename( {arr.dims[1] : 'lev'} )
    arr.rename( {arr.dims[2] : 'lat'} )

    if 'coordinates' in arr.attrs:
        arr.attrs.pop('coordinates')

    if arr.coords['lev'].attrs['units'] == 'centimeters' or arr.coords['lev'].attrs['units'] == 'cm':   # line 475 of ncl version
        lev = arr.coords['lev']
        lev.attrs['units'] = "m"
        lev = lev/100.
        ##RLB do thes next two statements carry over automatically?
        ## lev&lev = lev
        ## arr&lev = lev


    if arr.coords['lev'][2] < 0:   # check for negative levels     line 488 ncl script
        arr.coords['lev'] = arr.coords['lev'] * -1.
        if [True for l in arr.coords['lev'] if l < 0]:
            print("Error detected in MOC level sign conversion")
            print(arr.coords['lev'])
        arr.coords['lev'].attrs['positive'] = 'down'
        ##RLB again, are these necessary?
        ##lev&lev = lev
        ##arr&lev = lev

    if 'valid_range' in arr.attrs:    # check to make sure data is in valid range. Reset to stay within the valid range if needed.
        arr.data = xr.where(arr.data < arr.attrs['valid_range'][0], arr.attrs['valid_range'][0], arr.data)
        arr.data = xr.where(arr.data > arr.attrs['valid_range'][1], arr.attrs['valid_range'][1], arr.data)

    if [True for a in arr.data if abs(a) >= 1.e20]:   # check for inf values or values way out of range, reset to _FillValue.
        print(f'Values greater than 1.e20 or less than -1.e20 detected in {zpath}, resetting to _FillValue')
        arr.data = xr.where(abs(arr.data) >= 1.e20, arr.attrs['_FillValue'], arr.data)

    if yearS < int(sydata) or yearE > int(eydata):    # line 517 of ncl script
        print(f'Requested {yearS}-{yearE} time span is outside the input file {zpath} time span of {sydata}-{eydata}')
        print('Setting array to all missing')
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lev_lat")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        timeT =  ncl.yyyymm_time(sydata, eydata, "integer")
        time = timeT[slice(sydata*100+smdata, eydata*100+emdata)]
        if 'time' in arr.coords:
            arr.drop('time')
        dimz = arr.values.shape

        if dimz[0] == time.values.shape[0]:
            #### THIS MAY NOT BE RIGHT -- MAY NEED TO REASSIGN ARR RESULTS OF ASSIGN_COORDS
            arr.assign_coords(time=time)
        else:
            print('Possible mismatch detected between time specified in file name and file variables, setting array to missing')
            print(f'File = {zpath}')
            print(f'Read from file name: {min(time)}-{max(time)}')
            arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lev_lat")
            sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
            eydata = yearE     # as data array is totally missing..
            smdata = 1
            emdata = 12

    if arr.coords['lat'][0] >= 0:   # line 553 in ncl script
        farr = arr[slice(yearS*100+1, yearE*100+12),:,::-1]   # flip the latitudes
    else:
        farr = arr[slice(yearS*100+1, yearE*100+12),:,:]

    mocheck = np.array((yearS*100+1)-min(farr.coords['time'], (yearE*100+12) - max(farr.coords['time'])))
    if [True for mon in mocheck if mon != 0]:    # previously: if (mod(dimsizes(farr&time),12).ne.0) then
        if mocheck[0] != 0:
            print("First requested year is incomplete")
        if mocheck[1] != 0:
            print("Last requested year is incomplete")
        print(f'Incomplete data year(s) requested for file {zpath}, printing out time and creating blank array')
        print(f'Time requested: {yearS}-{yearE}')
        print(farr.coords['time'])
        farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    # check units for MOC array. CMIP5 = "kg s-1"  CMIP3 = "m3 s-1" CCSM3 = "Sverdrups" CCSM4 = "Sverdrups"

    if farr.attrs['units'] == 'Sverdrups':    # line 579 of ncl script
        farr.attrs['units'] = "Sv"
    if farr.attrs['units'] == "kg s-1" or farr.attrs['units'] == "KG S-1" or  \
            farr.attrs['units'] == "kg/s" or farr.attrs['units'] == "KG/S":     # 1 Sv = 1.e9 kg/s
        farr.values =  farr.values/1.e9      ###RLB -- original strips metadata???
        farr.attrs['units'] = "Sv"
    if farr.attrs['units'] == 'm3 s-1' or farr.attrs['units'] == 'M3 S-1' or \
            farr.attrs['units'] == 'm3/s' or farr.attrs['units'] == 'M3/S':     # 1 Sv = 1.e6 m3/s
        farr.values = farr.values/1.e6       ###RLB -- original strips metadata???
        farr.attrs['units']= "Sv"

    date = farr.coords['time']         # switch time to be CF-conforming
    yyyy = date / 100
    mm = date - (yyyy * 100)
    days = (ncl.days_in_month(yyyy, mm)) / 2
    hms = days
    hms = 0  # hours, minutes, seconds all the same (=0)
    time = ncl.cd_inv_calendar(yyyy, mm, days, hms, hms, hms, "months since " + min(yyyy) + "-01-15 00:00:00", 0)
    time.attrs['long_name'] = 'Time'
    time.attrs['standard_name'] = 'time'
    time.attrs['actual_range'] = np.array(min(time), max(time))
    time.rename({time.dims[0]: 'time'})
    time.coords['time'] = time
    farr.coords['time'] = time
    return farr

def data_read_in_ice(zpath, vn, yearS, yearE):
    '''
    read in ice data from given files
    assign time coordinate variables, check for issues with the array, assign _FillValue (if needed)
    assign dimension names (for ease-of-use), check and modify units
    '''

    arr = None

    if vn == 'aice_nh':            # line 620 of ncl script
        vname = ("aice_nh","aice","sic","SIC","CN","ice","icec")
    if vn == 'aice_sh':
        vname = ("aice_sh","aice","sic","SIC","CN","ice","icec")

    if zpath is None:
        print('File missing, creating blank array of data. View {vn} namelist for details.')
        arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        tfiles = glob.glob(zpath)           # line 638 of ncl script
        c = xr.open_mfdataset(tfiles)
        for v in vname:
            if v in c:
                if v.name == 'aice_nh' or v.name == 'aice_sh' or v.name == 'aice':   # CCSM/CESM file
                    arr = c[v]
                    if 'coordinates' in arr.attrs:
                        strarr = arr.attrs['coordinates'].split()
                        if 'TLON' in strarr:     # CESM longitude 2D coordinate
                            if len(c.TLON.dims) == 3:
                                arr.attrs['lon2d'] = c.TLON[0,:,:]
                        else:
                            arr.attrs['lon2d'] = c.TLON
                        if 'TLAT' in strarr:     # CESM latitude 2D coordinate
                            if len(c.TLAT.dims) == 3:
                                arr.attrs['lat2d'] = c.TLAT[0,:,:]
                            else:
                                arr.attrs['lat2d'] = c.TLAT

                    if 'cell_measures' in arr and 'tarea' in c:    # if an attribute named cell_measures exists, and tarea is on file(0)
                        if arr.attrs['cell_measures'] == 'area: tarea':
                            arr.attrs['area'] = c.tarea.values.astype(arr.dtype)      # in units of m^2
                else:
                    if vname == 'CN':     # GFDL file       line 674 in ncl script
                        arrT = c.variables[vname]
                        arr = arrT.sum(1)
                        arr = xr.where(arr > 1,1,arr)   # optional
                    else:
                        arr = c.variables[vname]
                    if "coordinates" in arr.attrs:
                        strarr = arr.attrs['coordinates'].split()
                        if 'lon' in strarr:      #IPCC longitude 2D coordinate
                            arr.attrs['lon2d'] = c.lon
                        if 'lat' in strarr:      # IPCC latitude 2D coordinate
                            arr.attrs['lat2d'] = c.lat
                        if 'longitude' in strarr:      # NSIDC longitude 2D coordinate
                            arr.attrs['lon2d'] = c.longitude
                        if 'latitude' in strarr:       # NSIDC latitude 2D coordinate
                            arr.attrs['lat2d'] = c.latitude
#                       else:
#                           print("2D coordinates for ice data are not detected")

                    dir_name = tfiles[0].split('/')   # line 700 ncl script
                    if len(dir_name) >= 8:
                        dir_name_new = "/".join(dir_name[:5]) + '/fx/areacello/' + dir_name[7] + '/r0i0p0/*.nc'
                        ufile = glob.glob(dir_name_new)
                    else:
                        ufile = ''
                    if len(ufile) > 0:
                        d = xr.open_dataset(ufile)
                        arr.attrs['area'] = d.areacello.values.astype(arr.dtype.name)
                        dimQ = arr.shape
                        if arr.attrs['area'].size != (dimQ[1]*dimQ[2]):   # the dimension sizes of areacello
                            arr.attrs.pop('area')                         # do not match sizes of area j,i dimensions

                    if 'AREA' in c:    # check to see if there is an AREA array present and if so use it
                        areaT = c.AREA
                        if areaT.attr['units'] == 'km^2':
                            area_unit_km2_to_m2 = True
                            areaT = areaT*1000000.
                            areaT.attr['units'] = "m^2"
                        areaT.attr['_FillValue'] = 1.e20
                        arr.attr['area'] = areaT.astype(arr.values.dtype.name)
                        if "pole_hole_area" in areaT.attrs:    # format of ystart, yend, value, ystart, yend, value
                            try:
                                "area_unit_km2_to_m2"    # ensure this variable exists...
                                extra_area = areaT.attrs['pole_hole_area'].astype("float32")
                                extra_area[2::3] = extra_area[2::3]*1000000.   # convert pole hole area from km^2->m^2
                                arr.attrs['pole_hole_area'] = extra_area.astype(arr.values.dtype.name)
                            except NameError:
                                arr.attrs['pole_hole_area'] = areaT.attrs['pole_hole_area'].astype(arr.values.dtype.name)
                break

        cpathS = tfiles[0]    # line 745 of NCL script
        cpathE = tfiles[len(tfiles) - 1]
        ncharS = len(cpathS)
        ncharE = len(cpathE)
        sydata = cpathS[ncharS - 17:ncharS - 14]
        smdata = cpathS[ncharS - 13:ncharS - 12]
        eydata = cpathE[ncharE - 10:ncharE - 7]
        emdata = cpathE[ncharE - 6:ncharE - 5]

    if arr is None:
        print(f'Variable ({vn}) not found. Examine input file {zpath}. Creating empty array and continuing')
        arr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    if 'area' not in arr.attrs:   # calculate grid cell areas manually (not implemented)
        # print("Grid cell areas not found.")
        pass


    if arr.dtype == 'int16':       # line 765 of NCL script
     arr = arr.values.astype('float32')

    if "_FillValue" not in arr.attrs:     # assign _FillValue if one is not present
        if 'missing_value' in arr.attrs:
            arr.attrs['_FillValue'] = arr.attrs['missing_value']
        else:
            arr.attrs['_FillValue'] = ncl.default_fillvalues[arr.dtype.name]


    arr.rename({arr.dims[0] : 'time'})    # line 778 of NCL script
    arr.rename({arr.dims[1] : 'j'})
    arr.rename({arr.dims[2] : 'i'})

    if 'lat2d' not in arr.attrs:        # if latitudes are 1D, make sure latitudes run from south to north +
        if arr.coords['j'][0] >= 0:     # calculate area of 1D lat/lon arrays
            tarr = arr[:,::-1,:]
            arr = tarr

        if arr.coords['i'].min() >= 0 and arr.coords['i'].max() <= 360:
            fctr = 111120   # how many meters per degree of latitude (approximate)
            pi=4.*math.atan(1.0)
            rad=(pi/180.)
            lat = arr.coords['j'].astype('float32')
            dimlat = lat.sizes
            latr = np.array(dimlat,lat.dtype)

            for gg in range(0,dimlat-1):      # line 797 in NCL script
                if gg == 0:
                    latr[gg] = abs(-90-(lat[1]+lat[0])/2.)
                elif gg == (dimlat-1):
                    latr[gg] = abs(90 - (lat[dimlat-2]+lat[dimlat-1])/2.)
                else:
                    latr[gg]= abs((lat[gg-1]+lat[gg])/2. - (lat[gg]+lat[gg+1])/2.)

            lon = arr.coords['i'].astype('float32')
            dimlon = lon.shape
            lonr = np.array(dimlon,'float32')
            for gg in range(0,dimlon-1):
                if gg == 0:
                    lonr[gg] = abs( (lon[1]+lon[0])/2. - (((lon[dimlon-1]+(lon[0]+360))/2.)-360) )
                elif gg == (dimlon-1):
                    lonr[gg] = abs(((lon[dimlon-1]+(lon[0]+360))/2.) - (lon[gg-1]+lon[gg])/2.)
                else:
                    lonr[gg] = abs((lon[gg]+lon[gg+1])/2. - (lon[gg-1]+lon[gg])/2.)

            area = arr[0,:,:].astype('float32')
            area = area.attrs['_FillValue']
            area.attrs['long_name'] = "Area of grid box"
            area.attrs['units'] = "m2"

            #   printVarSummary(area)
            for ff in range(0,dimlat-1):          # line 828 in NCL script
                for gg in range(0,dimlon-1):
                    ###RLB the following line had (/ ...  /) around the expression ??
                    area[ff,gg] = (fctr*latr[ff])*(math.cos(rad*lat[ff])*lonr[gg]*fctr)     # cosine weighting
            # print("Total area = "+sum(area))
            arr.attrs['area'] = area.astype(arr.dtype.name)

    if 'is_all_missing' not in arr.attrs:    # erase data in hemisphere not specified via vn
        if 'lat2d' in arr.attr:
            tlat2 = ncl.conform(arr,arr.coords['lat2d'],(1,2))
            tlon2 = ncl.conform(arr,arr.coords['lon2d'],(1,2))
            if vn == 'aice_nh':
                arr = xr.where(tlat2>= 0, arr, arr.attrs['_FillValue'])
            if vn == 'aice_sh':
                arr = xr.where(tlat2 < 0, arr, arr.attrs['_FillValue'])
        else:
            if vn == 'aice_nh':
                arr.loc[:,:-1.,:] = arr.attrs['_FillValue']   ##RLB not so sure about coordinate indexing here?
            if vn == 'aice_sh':
                arr.loc[:,0:,:] = arr.attrs['_FillValue']     ##RLB not so sure about coordinate indexing here?

    if yearS < int(sydata) or yearE > int(eydata):    # line 860 of NCL script...
        print(f'Requested {yearS}-{yearE} time span is outside the input file {zpath} time span of {sydata}-{eydata}')
        print("Setting array to all missing")
        arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
        sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
        eydata = yearE     # as data array is totally missing..
        smdata = 1
        emdata = 12
    else:
        timeT =  ncl.yyyymm_time(sydata, eydata, "integer")
        time = timeT.loc[sydata*100+smdata, eydata*100+emdata]
        if 'time' in arr.coords:
            arr.coords.pop('time')
        dimz = arr.shape
        if dimz[0] == time.sizes:
            arr.coords['time'] = time
        else:
            print('Possible mismatch detected between time specified in file name and file variables, setting array to missing')
            print(f'File = {zpath}')
            print("Read from file name: "+min(time)+"-"+max(time))
            arr = create_empty_array(yearS,yearE,smdata,emdata,"time_lat_lon")
            sydata = yearS     # assign these variables based on yearS/yearE provided in namelist. Doesn't matter
            eydata = yearE     # as data array is totally missing..
            smdata = 1
            emdata = 12

    farr = arr.loc[yearS*100+1:yearE*100+12,:,:]
    #   printVarSummary(farr)

    mocheck = np.array((yearS*100+1)-min(farr.coords['time']),(yearE*100+12) - max(farr.coords['time']))  # line 897 of NCL script
    if [True for m in mocheck if m != 0]:   # previously: if (mod(dimsizes(farr&time),12).ne.0) then
        if mocheck[0] != 0:
          print("First requested year is incomplete")
        if mocheck[1] != 0:
            print("Last requested year is incomplete")
        print(f'Incomplete data year(s) requested for file {zpath}, printing out time and creating blank array')
        print(f'Time requested: {yearS}-{yearE}')
        print(f'From file: Times present from {farr.coords["time"].min()}-{farr.coords["time"].max()})')
        farr = create_empty_array(yearS,yearE,1,12,"time_lat_lon")

    if farr.attrs['units'] == "0-1" or farr.attrs['units'] == "1":    # GFDL units, NSIDC units
     farr.values = farr.values * 100      ###RLB -- was (/ farr*100. /)   metadata stripped?
     farr.attrs['units'] = "%"

    date = farr.coords['time']         # switch time to be CF-conforming
    yyyy = date / 100
    mm = date - (yyyy*100)
    days = (ncl.days_in_month(yyyy,mm))/2
    hms  = days
    hms = 0   # hours, minutes, seconds all the same (=0)
    time = ncl.cd_inv_calendar(yyyy,mm,days,hms,hms,hms,"months since "+min(yyyy)+"-01-15 00:00:00",0)
    time.attrs['long_name'] = "Time"
    time.attrs['standard_name'] = "time"
    time.rename( {time.dims[0] : 'time'}  )
    time.coords['time'] = time
    farr.coords['time'] = time
    return farr


def y_axis_check(temparr, tempres):    # line 940 of NCL script
    '''
      alters the formatting of the Y-axis

      not currently used
    '''
    minval = temparr.min()
    maxval = temparr.max()
    if minval > -1 and minval < 0 and maxval < 1 and maxval > 0:
        tempres.attrs['tmYLFormat'] = "0@;*.2f"
    else:
        tempres.attrs['tmYLFormat'] = "0@*+^sg"
    return tempres


def check_custom_climo(mn, startyear, endyear, climo_startyear, climo_endyear):  # line 956 of NCL script
    '''
    Check that the user-specified climatological period is within the time range of the data
    '''

    for gg in range(0,startyear.sizes[0]-1):
        if climo_startyear>= 0:    # exact years specified for climatological period
            if climo_startyear >= startyear[gg] and climo_endyear<= endyear[gg]:
                pass
            else:
                print('check_custom_climo: Warning! Beginning and/or ending of climatological period is outside time range of data.')
                print(f'Dataset: {mn}, years = {startyear[gg]}:{endyear[gg]}, set climatological period = {climo_startyear}:{climo_endyear}')
                print('The diagnostics package will proceed, but one or more dataset(s) will not have the full climatological period removed and/or the package may fail with the following message: fatal:NclOneDValGetRangeIndex: start coordinate index out of range.')
        else:                              # relative years specified for climatological period
            if (endyear[gg]-startyear[gg]+1) < (climo_endyear-climo_startyear+1):
                print('check_custom_climo: Warning! Beginning and/or ending of climatological period is outside time range of data.')
                print(f'Dataset: {mn}, years = {startyear[gg]}:{endyear[gg]}, set climatological period = {(endyear(gg)+climo_startyear)}:{(endyear(gg)+climo_endyear)}')
                print('The diagnostics package will proceed, but one or more dataset(s) will not have the full climatological period removed and/or the package may fail with the following message: fatal:NclOneDValGetRangeIndex: start coordinate index out of range.')
            if climo_startyear.abs() >= (endyear[gg]-startyear+1):
                print(f'check_custom_climo: Warning! Dataset: {mn}, climatology start year {(endyear[gg]+climo_startyear)} is outside of analysis time period ({startyear[gg]}:{endyear[gg]}), exiting script.')
                quit()


def isfilepresent2(fdpres):      # line 989 of NCL script
    '''
     In version 6.2.1 the behavior of isfilepresent switched, where only files readable by NCL return
     True. Previously if a file (or directory) simply existed, isfilepresent would return True. A new
     function has been created in v6.2.1, fileexists, that acts like the previous version of isfilepresent
     did. To compensate for this, check the NCL version number, and use isfilepresent/fileexists when
     appropriate.

     RLB: This function would seem to be superfluous now -- we'll implement it in Python terms...
    '''
    return os.path.exists(fdpres)
#   nclver = stringtochar(get_ncl_version())
#
#   num0 = toint(tostring(nclver(0)))
#   num1 = toint(tostring(nclver(2)))
#   num2 = toint(tostring(nclver(4)))
#   if (num0.le.5) then
#       ra = isfilepresent(fdpres)
#   end if
#   if (num0.eq.6) then
#       if (num1.le.1) then
#           ra = isfilepresent(fdpres)
#      end if
#      if (num1.eq.2) then
#           if (num2.eq.0) then
#               ra = isfilepresent(fdpres)
#           else
#               ra = fileexists(fdpres)
#           end if
#       end if
#       if (num1.ge.3) then
#         ra = fileexists(fdpres)
#       end if
#   end if
#   if (num0.ge.7) then
#       ra = fileexists(fdpres)
#   end if
#   return(ra)
#   delete([/nclver,num0,num1,ra/])


def table_link_setup(ipath, iname, ltxt):    # line 1024 of NCL script
    '''
    ; image name, along with link text
    '''
    if isfilepresent2(ipath+iname):
        otxt = f'a href="{iname}">{ltxt}</a>'
    else:
        otxt = ltxt

    return otxt


def gsn_panel2(wksp, plotp, lpl, panelres):   # line 1039 of NCL script
    '''
    ; checks to make sure at least one image is present in plot before paneling,
    ; thereby eliminating error message:
    ; Error: gsn_panel: all of the plots passed to gsn_panel appear to be invalid
    '''

    if [True for p in plotp if not None]:
        ncl.gsn_panel(wksp,plotp,lpl,panelres)


def eofunc_north2(eval, N, prinfo):     # line 1052 of NCL script
    '''
    ;
    ; North, G.R. et al (1982): Sampling Errors in the Estimation of Empirical Orthogonal Functions.
    ; Mon. Wea. Rev., 110, 699–706.
    ; doi: http://dx.doi.org/10.1175/1520-0493(1982)110<0699:SEITEO>2.0.CO;2
    ;
    ; Usage after 'eofunc'. Here ntim was used,
    ;             prinfo = True
    ;             sig    = eval_north(eof@eval, ntim, prinfo)
    ;
    ; Copied directly from v6.3.0 contributed.ncl for use in the package regardless of NCL version.
    ;
    '''

    neval   = eval.shape
    if neval[0] == 1:
        print("eofunc_north: neval=1, no testing can be performed")
        sig = xr.DataArray([True])
        sig.attrs['long_name'] = "EOF separation is not testable N=1"
        sig.attrs['N']         =  N
        return sig

    dlam    = eval * math.sqrt(2.0/N)   # eq 24
    low     = eval-dlam
    high    = eval+dlam

    sig     = xr.DataArray(np.array((eval.shape[0]), dtype=np.bool))
    sig = False  # default is not significantly separated

    # first and last eigenvalues are special cales

    if eval[0] > high[1]:
        sig[0] = True

    if eval[neval-1] < low[neval-2]:
        sig[neval-1] = True

    # loop over other eignvalues

    if N > 2:
        for n in range(1,neval-2):
            if eval[n] < low[n-1] and eval[n] > high[n+1]:
                sig[n] = True

    if prinfo:
        print(f'{dlam}   {low}   {eval}   {high}  {sig}')

    sig.attrs['long_name'] = "EOF separation"
    sig.attrs['N']         =  N
    return sig


def set_varAtts(zarr, loname, uni, com_cvdp):   # line 1116 of NCL script
    '''
    ; Standardize and set attributes for array. Remove NCL-added and superfluous attributes,
    ; set missing_value equal to _FillValue, provide options to set long_name, units, and comment_cvdp
    ; attributes. For last three inputs "" means leave as set and "delete" means delete the attribute.
    ; This function will be used immediately prior to an array being written to a netCDF file.
    ;
    '''

    if 'anomaly_op_ncl' in zarr.attrs:
        zarr.attrs.pop('anomaly_op_ncl')

    if 'average_op_ncl' in zarr.attrs:
        zarr.attrs.pop('average_op_ncl')

    if 'cell_measures' in zarr.attrs:
        zarr.attrs.pop('cell_measures')

    if 'cell_methods' in zarr.attr:
        zarr.attrs.pop('cell_methods')

    if 'lonFlip' in zarr.attrs:
        zarr.attr.pop('lonFlip')

    if 'runave_op_ncl' in zarr.attrs:
        zarr.attrs.pop('runave_op_ncl')

    if 'stddev_op_ncl' in zarr.attrs:
        zarr.attrs.pop('stddev_op_ncl')

    if 'sum_op_ncl' in zarr.attrs:
        zarr.attrs.pop('sum_op_ncl')

    if 'time' in zarr.attrs:
        zarr.attrs.pop('time')

    if 'wgt_areaave_op_ncl' in zarr.attrs:
        zarr.attrs.pop('wgt_areaave_op_ncl')

    if '_FillValue' in zarr.attrs:        # set missing_value = _FillValue
        zarr.attrs['missing_value'] = zarr.attrs['_FillValue']

    if loname == 'delete':
        zarr.attrs.pop('long_name')
    elif loname != "":    # "" = leave as is
        zarr.attrs['long_name'] = loname

    if uni == 'delete':
        zarr.attrs.pop('units')
    elif uni != "":    # "" = leave as is
        zarr.attrs['units'] = uni

    if com_cvdp == 'delete':
        zarr.attrs.pop('comment_cvdp')
    elif com_cvdp != "":   # "" = leave as is
        zarr.attrs['comment_cvdp'] = com_cvdp

    return zarr

