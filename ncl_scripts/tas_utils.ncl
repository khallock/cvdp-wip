;=================================================================================================
; Open an intermediate NetCDF file for writing data for plots.
;
undef("open_file_for_plot_data")
function open_file_for_plot_data(prefixName,scriptName)
local fileAtts
begin
  fileName = getenv("OUTDIR") + prefixName + ".nc"
  if(fileexists(fileName))
    system("rm -f " + fileName)
  end if
  f = addfile(fileName,"c")
  fileAtts               = True            ; assign file attributes
  fileAtts@title         = "Intermediate file created for " + scriptName
  fileAtts@Conventions   = "None"
  fileAtts@creation_date = systemfunc ("date")

  fileattdef(f,fileAtts )                  ; copy file attributes
  return(f)
end
