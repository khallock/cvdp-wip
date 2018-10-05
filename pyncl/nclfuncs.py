import xarray as xr
import numpy as np

# Maps numpy dtypes onto NCL's corresponding default-values...
default_fillvalues = {
    'float32' : 9.96921e+36,
    'float64' : 9.969209968386869e+36,
    'int16'   : -32767,
    'int32'   : -2147483647,
    'int64'   : -9223372036854775806,
    'string'  : 'missing',
}

def yyyymm_time(yrStrt, yrLast, t=int):
    ''' usage
          yyyymm = yyyymm_time(1800,2001,int|float)'''
    nmos = 12
    mons = np.arange(1,13)
    nyrs = yrLast-yrStrt+1
    ntim = nmos*nyrs
    timeType = int
    if t == float:
        timeType = float
    timeValsNP = np.empty(ntim, dtype=timeType)

    n = 0
    for yr in range(yrStrt, yrLast+1):
        timeValsNP[n:n+nmos] = yr*100 + mons
        n = n+nmos

    timeVals = xr.DataArray(timeValsNP, dims=('time'), coords={'time': timeValsNP},
                            attrs={'long_name' : 'time', 'units' : 'YYYYMM'})
    return timeVals

def lonFlip(arr):
    print('************************************')
    print('lonFlip is stubbed and unimplemented')
    print('************************************')

def days_in_month(arr1, arr2):
    print('******************************************')
    print('days_in_month is stubbed and unimplemented')
    print('******************************************')

def cd_inv_calendar(d1, d2, d3, d4, d5, d6, d7, d8):
    print('********************************************')
    print('cd_inv_calendar is stubbed and unimplemented')
    print('*******************************************')

def conform(arr1, arr2, arr3):
    print('****************************************************************************')
    print('conform is stubbed and unimplemented or need mapping to xarray functionality')
    print('****************************************************************************')

def gsn_panel(arr1, arr2, arr3):
    print('*****************************************************************')
    print('gsn_panel is stubbed and should be implemented in terms of PyNGL ')
    print('*****************************************************************')
