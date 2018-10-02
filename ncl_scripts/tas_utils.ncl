function open_file(prefixName,scriptName)
local fileAtts
begin
  fileName = getenv("OUTDIR") + prefixName + ".nc"
  if(fileexists(fileName))
;   print("rm -f " + fileName)
    system("rm -f " + fileName)
  end if
; print("Creating file '" + fileName+"'")
  f = addfile(fileName,"c")
  fileAtts               = True            ; assign file attributes
  fileAtts@title         = "Intermediate file created for " + scriptName
  fileAtts@Conventions   = "None"
  fileAtts@creation_date = systemfunc ("date")

  fileattdef(f,fileAtts )                  ; copy file attributes
  return(f)
end
