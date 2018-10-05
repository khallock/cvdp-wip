import pyncl.nclfuncs as ncl

#
# tests of yyyymm_time
#
t = ncl.yyyymm_time(1800,1900)
print(len(t))
print(t)

t = ncl.yyyymm_time(1800,1900, float)
print(len(t))
print(t.sel(time=slice(180001,180012)))
