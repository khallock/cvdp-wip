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

;=================================================================================================
; Create a contour plot
;
undef("create_contour_map_plot")
function create_contour_map_plot(wks[1]:graphic,data[*][*]:numeric,name[1]:string,\
         start_year,end_year,colormap_type)
local res
begin
     res = True
     res@mpProjection = "WinkelTripel"
     res@mpGeophysicalLineColor = "gray42"
     if (NhlClassName(wks).eq."imageWorkstationClass")
        res@mpGeophysicalLineThicknessF = 2.  
     else
        res@mpGeophysicalLineThicknessF = 1.  
     end if     
     res@mpPerimOn    = False
     res@mpGridLatSpacingF =  90            ; change latitude  line spacing
     res@mpGridLonSpacingF = 180.           ; change longitude line spacing
     res@mpGridLineColor   = "transparent"  ; trick ncl into drawing perimeter
     res@mpGridAndLimbOn   = True           ; turn on lat/lon lines  
     res@mpFillOn = False
     res@mpCenterLonF = 210.
     res@mpOutlineOn = True  
     res@gsnDraw      = False
     res@gsnFrame     = False
  
     res@cnLevelSelectionMode = "ExplicitLevels"
     if (colormap_type.eq.0) then
        res@cnLevels = (/-8,-6,-5,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,5,6,8/)
     end if
     if (colormap_type.eq.1) then
        res@cnLevels = (/-6,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,6/)
     end if

     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014
     res@gsnLeftString = start_year+"-"+end_year

     res@gsnRightString = data@units
     res@gsnCenterString = name

     plot = gsn_csm_contour_map(wks,data,res)
     return(plot)
end