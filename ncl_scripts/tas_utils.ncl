;=================================================================================================
; Open an intermediate NetCDF file for writing data for plots.
;
undef("open_file_for_plot_data")
function open_file_for_plot_data(prefixName,scriptName)
local fileAtts
begin
  setfileoption("nc","Format","NetCDF4")       ; Use NetCDF4 so we can write strings
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
function create_contour_map_plot(wks[1]:graphic,data[*][*]:numeric,names[*]:string,\
         syear[*],eyear[*],ee[1]:integer)
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
        res@cnFillPalette = "ncl_default"
        res@cnLevels = (/-8,-6,-5,-4,-3,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,3,4,5,6,8/)
     end if
     if (COLORMAP.eq.1) then
        res@cnFillPalette = "BlueDarkRed18"
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
     res@gsnLeftString = syear(ee)+"-"+eyear(ee)

     res@gsnRightString = data@units
     res@gsnCenterString = names(ee)

     plot = gsn_csm_contour_map(wks,data,res)
     return(plot)
end

;=================================================================================================
; Create an XY plot
;
undef("create_xy_plot")
function create_xy_plot(wks[1]:graphic,data:numeric,data_obs:numeric,units[1]:string,\
                        names[*]:string,syear[*],eyear[*],ee[1]:integer,plot_type)
local nsim, nyr, nyr_max, xyres, xyres2, SCALE_TIMESERIES, tttt, plot_overlay,
      obs_file_present
begin
     nsim = dimsizes(syear)
     nyr = eyear-syear+1
     nyr_max = max(nyr)
     obs_file_present = isfilepresent2("obs_trefht")

     SCALE_TIMESERIES = getenv_logical("SCALE_TIMESERIES")

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
     if (obs_file_present.and.ee.eq.0) then
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
     if (SCALE_TIMESERIES) then
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
     if (obs_file_present.and.ee.ge.1) then
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
     if (obs_file_present.and.ee.ge.1) then
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
  SCALE_TIMESERIES = getenv_logical("SCALE_TIMESERIES")

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
  if (SCALE_TIMESERIES) then
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

;=================================================================================================
; Plots data from tas.trends_timeseries.ncl
;
undef("plot_tas_trends_timeseries")
procedure plot_tas_trends_timeseries()
local f, syear, eyear, names, f, nsim, nyr, obs_file_present
begin
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  OUTDIR           = getenv("OUTDIR")
  obs_file_present = isfilepresent2("obs_trefht")

;---Open one of the NetCDF files to get the years
  f = addfile(OUTDIR+"tas.trends.seas.ee0.nc","r")
  syear = f->syear
  eyear = f->eyear
  names = f->names
  delete(f)

  nsim = dimsizes(syear)
  nyr = eyear-syear+1

  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_trends_djf = gsn_open_wks(wks_type,OUTDIR+"tas.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,OUTDIR+"tas.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,OUTDIR+"tas.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,OUTDIR+"tas.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,OUTDIR+"tas.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,OUTDIR+"tas.trends.mon")
    
  wks_aa_djf = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.djf")
  wks_aa_mam = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.mam")
  wks_aa_jja = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.jja")
  wks_aa_son = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.son")
  wks_aa_ann = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.ann")
  wks_aa_mon = gsn_open_wks(wks_type,OUTDIR+"tas.timeseries.mon")
  
  wks_rt_mon = gsn_open_wks(wks_type,OUTDIR+"tas.runtrend.mon")

  map_djf = new(nsim,"graphic")  
  map_mam = new(nsim,"graphic")  
  map_jja = new(nsim,"graphic")  
  map_son = new(nsim,"graphic")  
  map_ann = new(nsim,"graphic")  
  map_mon = new(nsim,"graphic")  
  xy_djf = new(nsim,"graphic")  
  xy_mam = new(nsim,"graphic")  
  xy_jja = new(nsim,"graphic")  
  xy_son = new(nsim,"graphic")  
  xy_ann = new(nsim,"graphic")  
  xy_mon = new(nsim,"graphic")  
  xy_rt_mon = new((/5,nsim/),"graphic")  
  
  if (obs_file_present) then
     xy_obs_djf = new(nsim,"graphic")  
     xy_obs_mam = new(nsim,"graphic")  
     xy_obs_jja = new(nsim,"graphic")  
     xy_obs_son = new(nsim,"graphic")  
     xy_obs_ann = new(nsim,"graphic")  
     xy_obs_mon = new(nsim,"graphic")  
     ncdf_aa_seas_obs = addfile(OUTDIR+"tas.aa.seas.obs.nc","r")
     ncdf_aa_ann_obs  = addfile(OUTDIR+"tas.aa.ann.obs.nc","r")
     ncdf_aa_mon_obs  = addfile(OUTDIR+"tas.aa.mon.obs.nc","r")

     tas_aa_seas_obs = ncdf_aa_seas_obs->tas_aa_seas_obs
     tas_aa_ann_obs = ncdf_aa_ann_obs->tas_aa_ann_obs
     tas_aa_mon_obs = ncdf_aa_mon_obs->tas_aa_mon_obs
  end if

  do ee = 0,nsim-1
     ncdf_trends_seas = addfile(OUTDIR+"tas.trends.seas.ee"+ee+".nc","r")
     ncdf_trends_ann  = addfile(OUTDIR+"tas.trends.ann.ee"+ee+".nc","r")
     ncdf_trends_mon  = addfile(OUTDIR+"tas.trends.mon.ee"+ee+".nc","r")
     ncdf_aa_seas     = addfile(OUTDIR+"tas.aa.seas.ee"+ee+".nc","r")
     ncdf_aa_ann      = addfile(OUTDIR+"tas.aa.ann.ee"+ee+".nc","r")
     ncdf_aa_mon      = addfile(OUTDIR+"tas.aa.mon.ee"+ee+".nc","r")
     ncdf_rt_mon      = addfile(OUTDIR+"tas.rt.mon.ee"+ee+".nc","r")

     tas_trends_seas := ncdf_trends_seas->tas_trends_seas
     tas_trends_ann := ncdf_trends_ann->tas_trends_ann
     tas_trends_mon := ncdf_trends_mon->tas_trends_mon
     tas_aa_seas := ncdf_aa_seas->tas_aa_seas
     tas_aa_ann := ncdf_aa_ann->tas_aa_ann
     tas_aa_mon := ncdf_aa_mon->tas_aa_mon
     tas_rt_mon := ncdf_rt_mon->tas_rt_mon

     map_djf(ee) = create_contour_map_plot(wks_trends_djf,tas_trends_seas(0,:,:),names,syear,eyear,ee)
     map_mam(ee) = create_contour_map_plot(wks_trends_mam,tas_trends_seas(1,:,:),names,syear,eyear,ee)
     map_jja(ee) = create_contour_map_plot(wks_trends_jja,tas_trends_seas(2,:,:),names,syear,eyear,ee)
     map_son(ee) = create_contour_map_plot(wks_trends_son,tas_trends_seas(3,:,:),names,syear,eyear,ee)
     map_ann(ee) = create_contour_map_plot(wks_trends_ann,tas_trends_ann,names,syear,eyear,ee)
     map_mon(ee) = create_contour_map_plot(wks_trends_mon,tas_trends_mon,names,syear,eyear,ee)

     xy_djf(ee) = create_xy_plot(wks_aa_djf,tas_aa_seas(0,:),tas_aa_seas_obs(0,:),tas_trends_seas@units,names,syear,eyear,ee,"seas")
     xy_mam(ee) = create_xy_plot(wks_aa_mam,tas_aa_seas(1,:),tas_aa_seas_obs(1,:),tas_trends_seas@units,names,syear,eyear,ee,"seas")
     xy_jja(ee) = create_xy_plot(wks_aa_jja,tas_aa_seas(2,:),tas_aa_seas_obs(2,:),tas_trends_seas@units,names,syear,eyear,ee,"seas")
     xy_son(ee) = create_xy_plot(wks_aa_son,tas_aa_seas(3,:),tas_aa_seas_obs(3,:),tas_trends_seas@units,names,syear,eyear,ee,"seas")
     xy_ann(ee) = create_xy_plot(wks_aa_ann,tas_aa_ann,tas_aa_ann_obs,tas_trends_ann@units,names,syear,eyear,ee,"ann")
     xy_mon(ee) = create_xy_plot(wks_aa_mon,tas_aa_mon,tas_aa_mon_obs,tas_trends_mon@units,names,syear,eyear,ee,"mon")
     
     do ff = 0,4
        if (.not.all(ismissing(tas_rt_mon(ff,:))))     
           xy_rt_mon(ff,ee) = create_xy_plot(wks_rt_mon,tas_rt_mon(ff,:),tas_aa_mon,tas_rt_mon@units,names,syear,eyear,ee,"rt")
        end if
     end do   
  end do
  
  create_trends_panel_plot(wks_trends_djf,map_djf,"TAS Trends (DJF)",nyr)
  delete(wks_trends_djf)
  
  create_trends_panel_plot(wks_trends_mam,map_mam,"TAS Trends (MAM)",nyr)
  delete(wks_trends_mam)
  
  create_trends_panel_plot(wks_trends_jja,map_jja,"TAS Trends (JJA)",nyr)
  delete(wks_trends_jja)

  create_trends_panel_plot(wks_trends_son,map_son,"TAS Trends (SON)",nyr)
  delete(wks_trends_son)

  create_trends_panel_plot(wks_trends_ann,map_ann,"TAS Trends (Annual)",nyr)
  delete(wks_trends_ann)

  create_trends_panel_plot(wks_trends_mon,map_mon,"TAS Trends (Monthly)",nyr)
  delete(wks_trends_mon)

  create_global_avg_panel_plot(wks_aa_djf,xy_djf,"TAS Global Average (DJF)",nyr)
  delete(wks_aa_djf)

  create_global_avg_panel_plot(wks_aa_mam,xy_mam,"TAS Global Average (MAM)",nyr)
  delete(wks_aa_mam)

  create_global_avg_panel_plot(wks_aa_jja,xy_jja,"TAS Global Average (JJA)",nyr)
  delete(wks_aa_jja)

  create_global_avg_panel_plot(wks_aa_son,xy_son,"TAS Global Average (SON)",nyr)
  delete(wks_aa_son)

  create_global_avg_panel_plot(wks_aa_ann,xy_ann,"TAS Global Average (Annual)",nyr)
  delete(wks_aa_ann)

  create_global_avg_panel_plot(wks_aa_mon,xy_mon,"TAS Global Average (Monthly)",nyr)
  delete(wks_aa_mon)

  do ff = 0,dimsizes(tas_rt_mon&rt_year)-1
    create_global_avg_panel_plot(wks_rt_mon,xy_rt_mon(ff,:),"TAS Running " + tas_rt_mon&rt_year(ff) + "yr Trend (Monthly)",nyr)
  end do
  delete(wks_rt_mon)
  
  delete([/map_djf,map_mam,map_jja,map_son,map_ann,map_mon,xy_djf,xy_mam,xy_jja,xy_son,xy_ann,xy_mon,xy_rt_mon/])
  if (obs_file_present) then
     delete([/xy_obs_mam,xy_obs_jja,xy_obs_son,xy_obs_ann,xy_obs_mon/])
  end if
  if (wks_type.eq."png") then  
     do gg = 1,5
        if (isfilepresent2(OUTDIR+"tas.runtrend.mon.00000"+gg+".png")) then
           system("mv "+OUTDIR+"tas.runtrend.mon.00000"+gg+".png "+OUTDIR+"tas."+tas_rt_mon&rt_year(gg-1)+"yr_runtrend.mon.png")
        end if
     end do
  else
     if (isfilepresent2(OUTDIR+"tas.runtrend.mon.ps")) then
        system("psplit "+OUTDIR+"tas.runtrend.mon.ps "+OUTDIR+"tas_rt")
        do gg = 1,5
           if (isfilepresent2(OUTDIR+"tas_rt000"+gg+".ps")) then
              system("mv "+OUTDIR+"tas_rt000"+gg+".ps "+OUTDIR+"tas."+tas_rt_mon&rt_year(gg-1)+"yr_runtrend.mon.ps")
           end if
        end do
        system("rm "+OUTDIR+"tas.runtrend.mon.ps")
     end if
  end if
end
