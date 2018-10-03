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
  fileAtts@title         = "Intermediate file created by " + scriptName
  fileAtts@Conventions   = "None"
  fileAtts@creation_date = systemfunc ("date")

  fileattdef(f,fileAtts )                  ; copy file attributes
  return(f)
end

;=================================================================================================
; Create a contour/map plot
;
undef("create_contour_map_plot")
function create_contour_map_plot(wks[1]:graphic,data[*][*]:numeric,name[1]:string,\
         syear,eyear)
local res
begin
     COLORMAP = toint(getenv("COLORMAP"))

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
     if (COLORMAP.eq.0) then
        res@cnLevels = (/-8,-6,-5,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,5,6,8/)
     end if
     if (COLORMAP.eq.1) then
        res@cnLevels = (/-6,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,6/)
     end if

     res@cnLineLabelsOn = False
     res@cnFillOn        = True
     res@cnLinesOn       = False
     res@lbLabelBarOn    = False

     if (COLORMAP.eq.0) then
       res@cnFillPalette = "ncl_default"
     end if
     if (COLORMAP.eq.1) then
       res@cnFillPalette = "BlueDarkRed18"
     end if

     res@gsnLeftStringOrthogonalPosF = -0.05
     res@gsnLeftStringParallelPosF = .005
     res@gsnRightStringOrthogonalPosF = -0.05
     res@gsnRightStringParallelPosF = 0.96
     res@gsnRightString = ""
     res@gsnLeftString = ""
     res@gsnLeftStringFontHeightF = 0.014
     res@gsnCenterStringFontHeightF = 0.018
     res@gsnRightStringFontHeightF = 0.014
     res@gsnLeftString = syear+"-"+eyear

     res@gsnRightString = data@units
     res@gsnCenterString = name

     plot = gsn_csm_contour_map(wks,data,res)
     return(plot)
end

;=================================================================================================
; Create an XY plot
;
undef("create_xy_plot")
function create_xy_plot(wks[1]:graphic,data:numeric,data_obs:numeric,units[1]:string,\
                        names[*]:string,syear,eyear,ee,plot_type)
local nsim, nyr, nyr_max, xyres, xyres2, SCALE_TIMESERIES, tttt, plot_overlay
begin
     nsim = dimsizes(syear)
     nyr = eyear-syear+1
     nyr_max = max(nyr)

     SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")

     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"

     if (NhlClassName(wks).eq."imageWorkstationClass")
        xyres@xyLineThicknessF = 4.
     else
        xyres@xyLineThicknessF = 2.
     end if
     if (isfilepresent2("obs_trefht").and.ee.eq.0) then
        xyres@xyLineColor = "black"
     else
        xyres@xyLineColor = "royalblue"
     end if
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.0125
        xyres@tmYLLabelFontHeightF = 0.0125
        xyres@gsnLeftStringFontHeightF = 0.017
        xyres@gsnRightStringFontHeightF = 0.013
     else
        xyres@tmXBLabelFontHeightF = 0.018
        xyres@tmYLLabelFontHeightF = 0.018
        xyres@gsnLeftStringFontHeightF = 0.024
        xyres@gsnRightStringFontHeightF = 0.020
     end if
     xyres@gsnLeftStringOrthogonalPosF = 0.025
     xyres@gsnRightStringOrthogonalPosF = xyres@gsnLeftStringOrthogonalPosF
     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.15
     if (SCALE_TIMESERIES.eq."True") then
        xyres@vpWidthF = 0.9*((nyr(ee)*1.)/nyr_max)
     else
        xyres@vpWidthF = 0.9
     end if
     xyres@gsnLeftString = ""
     xyres@gsnCenterString = ""
     xyres@gsnRightString = ""

     xyres@trXMinF = syear(ee)-.5
     if(plot_type.eq."mon".or.plot_type.eq."rt") then
       xyres@trXMaxF = eyear(ee)+1.5
     else
       xyres@trXMaxF = eyear(ee)+0.5
     end if
     xyres2 = xyres
     xyres2@xyLineColor = "gray60"
     xyres2@xyCurveDrawOrder = "PreDraw"

     xyres@gsnLeftString = names(ee)
     tttt = dtrend_msg(ispan(0,dimsizes(data&time)-1,1),data,False,True)
     if(plot_type.eq."rt")
       xyres@gsnRightString = units
     else if(plot_type.eq."mon")
       xyres@gsnRightString = decimalPlaces(tttt@slope*dimsizes(data&time),2,True)+units
     else
       xyres@gsnRightString = decimalPlaces(tttt@slope*nyr(ee),2,True)+units
     end if
     end if
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        xyres@trYMinF = min((/min(data),min(data_obs)/))-.01
        xyres@trYMaxF = max((/max(data),max(data_obs)/))+.01
     end if
     if(plot_type.eq."mon") then
       plot_base  = gsn_csm_xy(wks,fspan(syear(ee),eyear(ee)+.91667,dimsizes(data)),data,xyres)
     else if(plot_type.eq."rt") then
;
; Note: data_obs is not really observational data for the "rt" case. We needed a variable with
; a time coordinate array, so data_obs is being used as a placeholder for this.
;
       plot_base = gsn_csm_xy(wks,fspan(syear(ee),eyear(ee)+.91667,dimsizes(data_obs&time)),data,xyres)
     else
       plot_base = gsn_csm_xy(wks,ispan(syear(ee),eyear(ee),1),data,xyres)
     end if
     end if
     if(plot_type.eq."rt")
        return(plot_base)
     end if
     if (isfilepresent2("obs_trefht").and.ee.ge.1) then
        overlay_string = "plot_overlay_"+plot_type+"_"+ee
        if(plot_type.eq."mon") then
          plot_tmp = gsn_csm_xy(wks,fspan(data_obs@syear,data_obs@eyear+.91667,dimsizes(data_obs)),data_obs,xyres2)
        else
          plot_tmp = gsn_csm_xy(wks,ispan(data_obs@syear,data_obs@eyear,1),data_obs,xyres2)
        end if
        overlay(plot_base,plot_tmp)
        copy_VarAtts(plot_tmp,plot_base)     ; Make sure plot_tmp's graphics attributes live outside this function!
        plot_base@$overlay_string$ = plot_tmp
     end if
     return(plot_base)
end

;=================================================================================================
; Create a panel plot for trends
;
undef("create_trends_panel_plot")
procedure create_trends_panel_plot(wks,plots,title,nyr)
local panres, nrow, ncol, nsim
begin
  nsim = dimsizes(nyr)
  panres = True
  panres@gsnMaximize = True  
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelLabelBar = True
  panres@gsnPanelYWhiteSpacePercent = 3.0
  panres@pmLabelBarHeightF = 0.05
  panres@pmLabelBarWidthF = 0.65
  panres@lbTitleOn = False
  panres@lbBoxLineColor = "gray70"
  panres@lbLabelFontHeightF = 0.013
  if (nsim.le.4) then
     if (nsim.eq.1) then
        panres@gsnPanelMainFontHeightF = 0.022
        panres@gsnPanelBottom = 0.50
     else
        panres@gsnPanelMainFontHeightF = 0.0145
        panres@gsnPanelBottom = 0.50
     end if
  else
     panres@gsnPanelMainFontHeightF = 0.016
     panres@gsnPanelBottom = 0.05
  end if
  panres@lbLabelStride = 1
  panres@gsnPanelMainString = title
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel2(wks,plots,(/nrow,ncol/),panres)

end

;=================================================================================================
; Create a panel plot for global average
;
undef("create_global_avg_panel_plot")
procedure create_global_avg_panel_plot(wks,plots,title,nyr)
local panres, nsim, nrow, ncol, nyr_max, tt, lp, SCALE_TIMESERIES
begin
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")

  nsim = dimsizes(nyr)
  nyr_max = max(nyr)

  panres = True
  panres@gsnMaximize = True
  panres@gsnPaperOrientation = "portrait"
  panres@gsnPanelYWhiteSpacePercent = 3.0  
  if (nsim.le.5) then
     panres@gsnPanelMainFontHeightF = 0.024
  else
     panres@gsnPanelMainFontHeightF = 0.016
  end if
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres@gsnPanelScalePlotIndex = tt(0)
  end if
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)  
  end if
  panres@gsnPanelMainString = title
  gsn_panel2(wks,plots,lp,panres)  
end
