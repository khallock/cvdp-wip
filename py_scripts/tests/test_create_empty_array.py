import py_scripts.functions as nclfuncs

#
# create_empty_array()
#
b = nclfuncs.create_empty_array(1980, 1990, 3, 8, 'time_lat_lon')
print(b)
print(b.coords['time'])

b = nclfuncs.create_empty_array(1980, 1990, 3, 8, 'time_lev_lat')
print(b)
print(b.coords['lev'])

