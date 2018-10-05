# Concatenate all .nc files from same model/observational dataset
# into a single .nc file.

import os
import glob
import subprocess
import shutil

def run():

    print("Starting: ncfiles.append.ncl")

    OUTDIR = os.getenv("OUTDIR")
    o = os.getenv("OBS")

    if o is True:
        obsflag = True
    else:
        obsflag = False
    na = []

    f = open('namelist')       # line 22 NCL script
    for line in f:
        if line != "":
            na.append(line)
    f.close()
    nsims = len(na)

    names = []                  # line 36 NCL script
    syear = []
    eyear = []

    delim = '|'
    for gg in na:
        tokens = gg.split(delim)
        name = tokens[0].strip()
        name = name.replace(' ', '_').replace('/', '_').replace("'", '_').replace('(', '_').replace(')', '_')
        names.append(name)
        syear.append(tokens[2].strip())
        eyear.append(tokens[3].strip())

    for gg in range(len(na)):       # line 47 NCL script
        fils = glob.glob(OUTDIR + names[gg] + '*.nc')
        if len(fils) == 0:
            # print("NetCDF files not found for " + names + ", skipping appending")
            continue

        ofile = OUTDIR + names[gg] + '.cvdp_data.' + syear[gg] + '-' + eyear[gg] + '.nc'
        if len(fils) == 1:
            shutil.move(fils[0], ofile)
        else:
            # if file master is present append individual data files to file master.
            if os.path.exists(ofile):       # line 55 NCL script
                for filshh in fils:
                    if filshh != ofile:
                        subprocess.run(['ncks', '-A', '-h', filshh, ofile])

            else:
                # if file master is not present, append individual data files to last file in list,
                # and when done move the last file to be the master file
                for filshh in fils[:len(fils)]
                    subprocess.run(['ncks','-A', '-h', filshh, fils[len(fils)]] )
                shutil.move(fils[len(fils)], ofile)

            if len(fils[:len(fils)-1]) >= 2:
                ##RLB - am confused about what this is supposed to do?
                #os.remove()
                #if (dimsizes(fils(:dimf-2)).ge.2) then
                #   system("rm "+str_sub_str(str_join(fils(:dimf-2)," "),ofile,""))   ; remove each script's file, but do not remove the master file (if present)

        subprocess.run(['ncks', '-O', ofile, ofile])     # line 71 NCL script

    # ------------------------------------------------
    if obsflag:         # line 80 NCL script

        f = open('obs_maxnum')
        maxnumobs = int(f.read())
        f.close()

        namelist_files = ['psl', 'prect', 'trefht', 'ts', 'snowdp', 'moc', 'aice_nh', 'aice_sh']

        delim = "|"
        cntr = 0
        namesB = ['missing']
        namesB *= ((maxnumobs * len(namelist_files)) - 1)    # want this a specific length

        f = open('namelist')  # line 22 NCL script
        for line in f:
            if line != "":
             namesB.append(line)
        f.close()

     namesB = new(maxnumobs*dimsizes(namelist_files),string)
     do gg = 0,dimsizes(namelist_files)-1                    ; grab all observational dataset names from namelist_$var files
        na = asciiread("namelist_byvar/namelist_"+namelist_files(gg),(/maxnumobs/),"string")
        namesB(cntr:cntr+maxnumobs-1) = str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_sub_str(str_strip(str_get_field(na,1,delim))," ","_"),"/","_"),"'","_"),"(","_"),")","_")
        cntr = cntr+maxnumobs
        delete(na)
     end do

     namesB = where(namesB.eq."",namesB@_FillValue,namesB)     ; for blank names set them to _FillValue
     if (any(namesB.eq."missing")) then
        namesB(str_match_ind(namesB,"missing")) = namesB@_FillValue ; check for any names containing "missing", set to _FillValue
     end if
     delete([/delim,cntr,namelist_files/])

     do gg = 0,dimsizes(namesB)-1
        if (.not.ismissing(namesB(gg))) then
           fils = systemfunc("ls "+OUTDIR+namesB(gg)+".cvdp_data.*.nc 2> /dev/null")
           if (.not.ismissing(fils(0))) then
              dimf = dimsizes(fils)
              fil0 = tochar(fils(0))
              suffix = tostring(fil0(dimsizes(fil0)-12:dimsizes(fil0)-1))
              delete(fil0)
              ofi = OUTDIR+namesB(gg)+".cvdp_data."+suffix
              if (dimf.ge.2) then
                 if (isfilepresent(ofi)) then                   ; if file master is present append individual data files to file master.
                    do hh = 0,dimf-1
                       if (fils(hh).ne.ofi) then
                          system("ncks -A -h "+fils(hh)+" "+ofi)
                       end if
                    end do
                 else                                                        ; if file master is not present, append individual data files to last file in list,
                    do hh = 0,dimf-2                                         ; and when done move the last file to be the master file
                       system("ncks -A -h "+fils(hh)+" "+fils(dimf-1))
                    end do
                    system("mv "+fils(dimf-1)+" "+ofi)
                 end if

                 if (dimsizes(fils(:dimf-2)).ge.2) then
                    system("rm "+str_sub_str(str_join(fils(:dimf-2)," "),ofi,""))   ; remove each script's file, but do not remove the master file (if present)
                 end if
              else
                 if (fils(0).ne.ofi) then
                    system("mv "+fils(0)+" "+ofi)
                 end if
              end if
              system("ncks -O "+ofi+" "+ofi)   ; done to alphabetize output variable
              delete([/dimf,ofi/])
           else
;              print("NetCDF files not found for "+namesB(gg)+", skipping appending")
           end if
           delete(fils)
        end if
     end do
     delete([/namesB/])
  end if
  print("Finished: ncfiles.append.ncl")
end


if __name__ == '__main__':
    run()
