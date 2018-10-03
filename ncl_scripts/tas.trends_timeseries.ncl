; Calculates 2m air temperature global trends, running global trends and timeseries
;
; Variables used: tas
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$CVDP_SCRIPTS/functions.ncl"
load "$CVDP_SCRIPTS/tas_utils.ncl"

begin
  script_name = get_script_prefix_name()
  print("Starting: " + script_name)
  
  SCALE_TIMESERIES = getenv("SCALE_TIMESERIES")  
  OUTPUT_DATA      = getenv("OUTPUT_DATA")  
  PNG_SCALE        = tofloat(getenv("PNG_SCALE"))
  OPT_CLIMO        = getenv("OPT_CLIMO")
  CLIMO_SYEAR      = toint(getenv("CLIMO_SYEAR"))
  CLIMO_EYEAR      = toint(getenv("CLIMO_EYEAR"))
  OUTPUT_TYPE      = getenv("OUTPUT_TYPE") 
  COLORMAP         = getenv("COLORMAP")  
  OUTPUT_PLOT_DATA = getenv("OUTPUT_PLOT_DATA")  
  
  nsim = numAsciiRow("namelist_byvar/namelist_trefht")
  na = asciiread("namelist_byvar/namelist_trefht",(/nsim/),"string")
  names = new(nsim,"string")
  paths = new(nsim,"string")
  syear = new(nsim,"integer",-999)
  eyear = new(nsim,"integer",-999)
  delim = "|"

  do gg = 0,nsim-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
     paths(gg) = str_strip(str_get_field(na(gg),2,delim))
     syear(gg) = stringtointeger(str_strip(str_get_field(na(gg),3,delim)))
     eyear(gg) = stringtointeger(str_strip(str_get_field(na(gg),4,delim)))
  end do
  nyr = eyear-syear+1
  nyr_max = max(nyr)

  pi=4.*atan(1.0)
  rad=(pi/180.)
         
  wks_type = OUTPUT_TYPE
  if (wks_type.eq."png") then
     wks_type@wkWidth = 1500*PNG_SCALE
     wks_type@wkHeight = 1500*PNG_SCALE
  end if
  wks_trends_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.djf")
  wks_trends_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.mam")
  wks_trends_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.jja")
  wks_trends_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.son")
  wks_trends_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.ann")
  wks_trends_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.trends.mon")
    
  wks_aa_djf = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.djf")
  wks_aa_mam = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.mam")
  wks_aa_jja = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.jja")
  wks_aa_son = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.son")
  wks_aa_ann = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.ann")
  wks_aa_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.timeseries.mon")
  
  wks_rt_mon = gsn_open_wks(wks_type,getenv("OUTDIR")+"tas.runtrend.mon")

  if (COLORMAP.eq.0) then
     gsn_define_colormap(wks_trends_djf,"ncl_default")   
     gsn_define_colormap(wks_trends_mam,"ncl_default")  
     gsn_define_colormap(wks_trends_jja,"ncl_default") 
     gsn_define_colormap(wks_trends_son,"ncl_default") 
     gsn_define_colormap(wks_trends_ann,"ncl_default") 
     gsn_define_colormap(wks_trends_mon,"ncl_default") 
     gsn_define_colormap(wks_aa_djf,"ncl_default")   
     gsn_define_colormap(wks_aa_mam,"ncl_default")  
     gsn_define_colormap(wks_aa_jja,"ncl_default") 
     gsn_define_colormap(wks_aa_son,"ncl_default") 
     gsn_define_colormap(wks_aa_ann,"ncl_default") 
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_rt_mon,"ncl_default") 
  end if
  if (COLORMAP.eq.1) then
     gsn_define_colormap(wks_trends_djf,"BlueDarkRed18")     
     gsn_define_colormap(wks_trends_mam,"BlueDarkRed18")    
     gsn_define_colormap(wks_trends_jja,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_son,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_ann,"BlueDarkRed18")   
     gsn_define_colormap(wks_trends_mon,"BlueDarkRed18")   
     gsn_define_colormap(wks_aa_djf,"ncl_default")   
     gsn_define_colormap(wks_aa_mam,"ncl_default")  
     gsn_define_colormap(wks_aa_jja,"ncl_default") 
     gsn_define_colormap(wks_aa_son,"ncl_default") 
     gsn_define_colormap(wks_aa_ann,"ncl_default") 
     gsn_define_colormap(wks_aa_mon,"ncl_default") 
     gsn_define_colormap(wks_rt_mon,"ncl_default") 
  end if
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
  
  if (isfilepresent2("obs_trefht")) then
     xy_obs_djf = new(nsim,"graphic")  
     xy_obs_mam = new(nsim,"graphic")  
     xy_obs_jja = new(nsim,"graphic")  
     xy_obs_son = new(nsim,"graphic")  
     xy_obs_ann = new(nsim,"graphic")  
     xy_obs_mon = new(nsim,"graphic")  
     if (OUTPUT_PLOT_DATA.eq."True") then
       ncdf_aa_seas_obs = open_file_for_plot_data("tas.aa.seas.obs",script_name)
       ncdf_aa_ann_obs  = open_file_for_plot_data("tas.aa.ann.obs",script_name)
       ncdf_aa_mon_obs  = open_file_for_plot_data("tas.aa.mon.obs",script_name)
     end if
  end if
  do ee = 0,nsim-1
     if (OUTPUT_PLOT_DATA.eq."True") then
       ncdf_trends_seas = open_file_for_plot_data("tas.trends.seas.ee"+ee,script_name)
       ncdf_trends_ann  = open_file_for_plot_data("tas.trends.ann.ee"+ee,script_name)
       ncdf_trends_mon  = open_file_for_plot_data("tas.trends.mon.ee"+ee,script_name)
       ncdf_aa_seas     = open_file_for_plot_data("tas.aa.seas.ee"+ee,script_name)
       ncdf_aa_ann      = open_file_for_plot_data("tas.aa.ann.ee"+ee,script_name)
       ncdf_aa_mon      = open_file_for_plot_data("tas.aa.mon.ee"+ee,script_name)
       ncdf_rt_mon      = open_file_for_plot_data("tas.rt.mon.ee"+ee,script_name)
     end if

     tas = data_read_in(paths(ee),"TREFHT",syear(ee),eyear(ee))    ; read in data, orient lats/lons correctly, set time coordinate variable up
     if (isatt(tas,"is_all_missing")) then
        delete(tas)
        continue
     end if 
     if (OPT_CLIMO.eq."Full") then
        tas = rmMonAnnCycTLL(tas)
     else
        check_custom_climo(names(ee),syear(ee),eyear(ee),CLIMO_SYEAR,CLIMO_EYEAR)
        temp_arr = tas
        delete(temp_arr&time)
        temp_arr&time = cd_calendar(tas&time,-1)
        if (CLIMO_SYEAR.lt.0) then
           climo = clmMonTLL(temp_arr({(eyear(ee)+CLIMO_SYEAR)*100+1:(eyear(ee)+CLIMO_EYEAR)*100+12},:,:))       
        else
           climo = clmMonTLL(temp_arr({CLIMO_SYEAR*100+1:CLIMO_EYEAR*100+12},:,:))                 
        end if
        delete(temp_arr)
        tas   = calcMonAnomTLL(tas,climo) 
        delete(climo)
     end if
  
     coswgt=cos(rad*tas&lat)
     coswgt!0 = "lat"
     coswgt&lat= tas&lat
     
     tas_aa_mon = wgt_areaave_Wrap(tas,coswgt,1.0,0)
     tttt = dtrend_msg_n(ispan(0,dimsizes(tas&time)-1,1),tas,False,True,0)
     tas_trends_mon = tas(0,:,:)
     tas_trends_mon = (/ onedtond(tttt@slope, (/dimsizes(tas&lat),dimsizes(tas&lon)/) ) /)
     tas_trends_mon = tas_trends_mon*dimsizes(tas&time)
     tas_trends_mon@units = tas@units+" "+nyr(ee)+"yr~S~-1~N~"
     delete(tttt)
     
     tas_seas = runave_n_Wrap(tas,3,0,0)
     tas_seas(0,:,:) = (/ dim_avg_n(tas(:1,:,:),0) /)
     tas_seas(dimsizes(tas&time)-1,:,:) = (/ dim_avg_n(tas(dimsizes(tas&time)-2:,:,:),0) /)
     tas_ann = runave_n_Wrap(tas,12,0,0)
     delete(tas)
     
     tas_trends_seas = tas_seas(:3,:,:)
     tas_trends_seas = tas_trends_seas@_FillValue
     tas_trends_ann  = tas_trends_seas(0,:,:)
     tas_aa_seas = new((/4,nyr(ee)/),typeof(tas_seas))
     tas_aa_seas!0 = "season"
     tas_aa_seas!1 = "time"
     tas_aa_seas&time = ispan(syear(ee),eyear(ee),1)
     tas_aa_seas&time@units = "YYYY"
     tas_aa_seas&time@long_name = "time"
     tas_aa_ann = tas_aa_seas(0,:)
     do ff = 0,4
        if (ff.le.3) then
           tarr = tas_seas(ff*3::12,:,:)     
        end if  
        if (ff.eq.4) then
           tarr = tas_ann(5::12,:,:)
        end if
        tttt = dtrend_msg_n(ispan(0,dimsizes(tarr&time)-1,1),tarr,False,True,0)   
        if (ff.le.3) then
           tas_trends_seas(ff,:,:) = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           tas_aa_seas(ff,:) = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        if (ff.eq.4) then
           tas_trends_ann = (/ onedtond(tttt@slope, (/dimsizes(tarr&lat),dimsizes(tarr&lon)/) ) /)
           tas_aa_ann = (/ wgt_areaave(tarr,coswgt,1.0,0) /)
        end if
        delete([/tarr,tttt/])        
     end do
     tas_trends_seas = tas_trends_seas*nyr(ee)
     tas_trends_seas@units = tas_seas@units+" "+nyr(ee)+"yr~S~-1~N~"
     tas_trends_ann = tas_trends_ann*nyr(ee)
     tas_trends_ann@units = tas_ann@units+" "+nyr(ee)+"yr~S~-1~N~"         
     delete([/tas_seas,tas_ann,coswgt/])
     
     if (isfilepresent2("obs_trefht").and.ee.eq.0) then
        tas_aa_seas@syear = syear(ee)
        tas_aa_seas@eyear = eyear(ee)
        tas_aa_mon@syear = syear(ee)
        tas_aa_mon@eyear = eyear(ee)
        tas_aa_ann@syear = syear(ee)
        tas_aa_ann@eyear = eyear(ee)
        tas_aa_seas_obs = tas_aa_seas
        tas_aa_mon_obs  = tas_aa_mon
        tas_aa_ann_obs  = tas_aa_ann
     end if
     
     dimT = dimsizes(tas_aa_mon)      ; calculate running trends from the monthly data
     rt_nyr = (/8,10,12,14,16/)
     tas_rt_mon = new((/5,dimT/),typeof(tas_aa_mon))
     tas_rt_mon!0 = "rt_year"
     tas_rt_mon!1 = "time"
     tas_rt_mon&time = tas_aa_mon&time    
     tas_rt_mon&rt_year = rt_nyr
     copy_VarAtts(tas_aa_mon,tas_rt_mon)
     tas_rt_mon@long_name =  tas_rt_mon@long_name+" global average running trend"
     do ff = 0,dimsizes(rt_nyr)-1
        incr = rt_nyr(ff)*12
        do gg = 0,dimT-incr-1
           tttt = dtrend_msg(ispan(0,incr-1,1),tas_aa_mon(gg:gg+incr-1),False,True)
           tas_rt_mon(ff,gg) = (/ tttt@slope*incr /)
           delete(tttt) 
        end do
     end do
     delete([/dimT,incr/])  

     if (OUTPUT_DATA.eq."True") then
        modname = str_sub_str(names(ee)," ","_")
        bc = (/"/","'","(",")"/)
        do gg = 0,dimsizes(bc)-1
           modname = str_sub_str(modname,bc(gg),"_")
        end do
        fn = getenv("OUTDIR")+modname+".cvdp_data.tas.trends_timeseries."+syear(ee)+"-"+eyear(ee)+".nc"
        if (.not.isfilepresent2(fn)) then
           z = addfile(fn,"c")
           z@source = "NCAR Climate Analysis Section's Climate Variability Diagnostics Package v"+getenv("VERSION")
           z@notes = "Data from "+names(ee)+" from "+syear(ee)+"-"+eyear(ee)
           if (OPT_CLIMO.eq."Full") then
              z@climatology = syear(ee)+"-"+eyear(ee)+" climatology removed prior to all calculations (other than means)"
           else
              if (CLIMO_SYEAR.lt.0) then
                 z@climatology = (eyear(ee)+CLIMO_SYEAR)+"-"+(eyear(ee)+CLIMO_EYEAR)+" climatology removed prior to all calculations (other than means)"
              else
                 z@climatology = CLIMO_SYEAR+"-"+CLIMO_EYEAR+" climatology removed prior to all calculations (other than means)"
              end if
           end if
           z@Conventions = "CF-1.6" 
        else
           z = addfile(fn,"w")
        end if
        tas_aa_seas2 = tas_aa_seas
        tas_aa_seas2!1 = "TIME"
        tas_aa_seas2&TIME = ispan(syear(ee),eyear(ee),1)
        tas_aa_seas2&TIME@units = "YYYY"
        tas_aa_seas2&TIME@long_name = "time"
        tas_aa_ann2 = tas_aa_ann
        tas_aa_ann2!0 = "TIME"
        tas_aa_ann2&TIME = ispan(syear(ee),eyear(ee),1)
        tas_aa_ann2&TIME@units = "YYYY"
        tas_aa_ann2&TIME@long_name = "time"
        z->tas_global_avg_mon = set_varAtts(tas_aa_mon,"tas global area-average (monthly)","C","")
        z->tas_global_avg_djf = set_varAtts(tas_aa_seas2(0,:),"tas global area-average (DJF)","C","")
        z->tas_global_avg_mam = set_varAtts(tas_aa_seas2(1,:),"tas global area-average (MAM)","C","")
        z->tas_global_avg_jja = set_varAtts(tas_aa_seas2(2,:),"tas global area-average (JJA)","C","")
        z->tas_global_avg_son = set_varAtts(tas_aa_seas2(3,:),"tas global area-average (SON)","C","")
        z->tas_global_avg_ann = set_varAtts(tas_aa_ann2,"tas global area-average (annual)","C","")
        z->$("tas_global_avg_runtrend_"+rt_nyr(0)+"yr")$ = set_varAtts(tas_rt_mon(0,:),"tas global area-average "+rt_nyr(0)+"yr running trend","","") 
        z->$("tas_global_avg_runtrend_"+rt_nyr(1)+"yr")$ = set_varAtts(tas_rt_mon(1,:),"tas global area-average "+rt_nyr(1)+"yr running trend","","") 
        z->$("tas_global_avg_runtrend_"+rt_nyr(2)+"yr")$ = set_varAtts(tas_rt_mon(2,:),"tas global area-average "+rt_nyr(2)+"yr running trend","","") 
        z->$("tas_global_avg_runtrend_"+rt_nyr(3)+"yr")$ = set_varAtts(tas_rt_mon(3,:),"tas global area-average "+rt_nyr(3)+"yr running trend","","") 
        z->$("tas_global_avg_runtrend_"+rt_nyr(4)+"yr")$ = set_varAtts(tas_rt_mon(4,:),"tas global area-average "+rt_nyr(4)+"yr running trend","","") 
        z->tas_trends_djf     = set_varAtts(tas_trends_seas(0,:,:),"tas linear trends (DJF)","","")
        z->tas_trends_mam     = set_varAtts(tas_trends_seas(1,:,:),"tas linear trends (MAM)","","")
        z->tas_trends_jja     = set_varAtts(tas_trends_seas(2,:,:),"tas linear trends (JJA)","","")
        z->tas_trends_son     = set_varAtts(tas_trends_seas(3,:,:),"tas linear trends (SON)","","")
        z->tas_trends_ann     = set_varAtts(tas_trends_ann,"tas linear trends (annual)","","")
        z->tas_trends_mon     = set_varAtts(tas_trends_mon,"tas linear trends (monthly)","","")
        delete(z)
        delete([/tas_aa_seas2,tas_aa_ann2/])
     end if    
;========================================================================
     map_djf(ee) = create_contour_map_plot(wks_trends_djf,tas_trends_seas(0,:,:),names(ee),syear(ee),eyear(ee))
     map_mam(ee) = create_contour_map_plot(wks_trends_mam,tas_trends_seas(1,:,:),names(ee),syear(ee),eyear(ee))
     map_jja(ee) = create_contour_map_plot(wks_trends_jja,tas_trends_seas(2,:,:),names(ee),syear(ee),eyear(ee))
     map_son(ee) = create_contour_map_plot(wks_trends_son,tas_trends_seas(3,:,:),names(ee),syear(ee),eyear(ee))
     map_ann(ee) = create_contour_map_plot(wks_trends_ann,tas_trends_ann,names(ee),syear(ee),eyear(ee))
     map_mon(ee) = create_contour_map_plot(wks_trends_mon,tas_trends_mon,names(ee),syear(ee),eyear(ee))

;---Write intermediate data used for contour plots
     if (OUTPUT_PLOT_DATA.eq."True") then
       ncdf_trends_seas->tas_trends_seas = tas_trends_seas
       ncdf_trends_ann->tas_trends_ann   = tas_trends_ann
       ncdf_trends_mon->tas_trends_mon   = tas_trends_mon
     end if
     xy_djf(ee) = create_xy_plot(wks_aa_djf,tas_aa_seas(0,:),tas_aa_seas_obs(0,:),tas_trends_seas@units,names,syear,eyear,ee,nsim,"seas")
     xy_mam(ee) = create_xy_plot(wks_aa_mam,tas_aa_seas(1,:),tas_aa_seas_obs(1,:),tas_trends_seas@units,names,syear,eyear,ee,nsim,"seas")
     xy_jja(ee) = create_xy_plot(wks_aa_jja,tas_aa_seas(2,:),tas_aa_seas_obs(2,:),tas_trends_seas@units,names,syear,eyear,ee,nsim,"seas")
     xy_son(ee) = create_xy_plot(wks_aa_son,tas_aa_seas(3,:),tas_aa_seas_obs(3,:),tas_trends_seas@units,names,syear,eyear,ee,nsim,"seas")
     xy_ann(ee) = create_xy_plot(wks_aa_ann,tas_aa_ann,tas_aa_ann_obs,tas_trends_ann@units,names,syear,eyear,ee,nsim,"ann")
     xy_mon(ee) = create_xy_plot(wks_aa_mon,tas_aa_mon,tas_aa_mon_obs,tas_trends_mon@units,names,syear,eyear,ee,nsim,"mon")
     
     do ff = 0,4
        if (.not.all(ismissing(tas_rt_mon(ff,:))))     
           xy_rt_mon(ff,ee) = create_xy_plot(wks_rt_mon,tas_rt_mon(ff,:),tas_aa_mon,tas_rt_mon@units,names,syear,eyear,ee,nsim,"rt")
        end if
     end do   
         
;---Write intermediate data used for XY plots
     if (OUTPUT_PLOT_DATA.eq."True") then
       ncdf_aa_seas->tas_aa_seas = tas_aa_seas
       ncdf_aa_ann->tas_aa_ann   = tas_aa_ann
       ncdf_aa_mon->tas_aa_mon   = tas_aa_mon
       ncdf_rt_mon->tas_rt_mon   = tas_rt_mon
     end if
     delete([/tas_trends_seas,tas_trends_ann,tas_trends_mon/])  
     delete([/tas_aa_seas,tas_aa_mon,tas_aa_ann,tas_rt_mon/])
  end do
  if (isfilepresent2("obs_trefht")) then
;---Write intermediate obs data used for XY plots
     if (OUTPUT_PLOT_DATA.eq."True") then
       ncdf_aa_seas_obs->tas_aa_seas_obs = tas_aa_seas_obs
       ncdf_aa_ann_obs->tas_aa_ann_obs   = tas_aa_ann_obs
       ncdf_aa_mon_obs->tas_aa_mon_obs   = tas_aa_mon_obs
     end if
     delete([/tas_aa_seas_obs,tas_aa_mon_obs,tas_aa_ann_obs/])
  end if
  
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
        panres@txFontHeightF = 0.022
        panres@gsnPanelBottom = 0.50
     else
        panres@txFontHeightF = 0.0145
        panres@gsnPanelBottom = 0.50
     end if
  else
     panres@txFontHeightF = 0.016
     panres@gsnPanelBottom = 0.05
  end if
  panres@lbLabelStride = 1
  
  panres@txString = "TAS Trends (DJF)"
  ncol = floattointeger(sqrt(nsim))
  nrow = (nsim/ncol)+mod(nsim,ncol)  
  gsn_panel2(wks_trends_djf,map_djf,(/nrow,ncol/),panres)
  delete(wks_trends_djf)
  
  panres@txString = "TAS Trends (MAM)"
  gsn_panel2(wks_trends_mam,map_mam,(/nrow,ncol/),panres)
  delete(wks_trends_mam)
  
  panres@txString = "TAS Trends (JJA)"
  gsn_panel2(wks_trends_jja,map_jja,(/nrow,ncol/),panres)
  delete(wks_trends_jja)
  
  panres@txString = "TAS Trends (SON)"
  gsn_panel2(wks_trends_son,map_son,(/nrow,ncol/),panres)
  delete(wks_trends_son)
  
  panres@txString = "TAS Trends (Annual)"
  gsn_panel2(wks_trends_ann,map_ann,(/nrow,ncol/),panres)
  delete(wks_trends_ann)
  
  panres@txString = "TAS Trends (Monthly)"
  gsn_panel2(wks_trends_mon,map_mon,(/nrow,ncol/),panres)
  delete(wks_trends_mon)
  delete(panres)
  
  panres2 = True
  panres2@gsnMaximize = True
  panres2@gsnPaperOrientation = "portrait"
  panres2@gsnPanelYWhiteSpacePercent = 3.0  
  if (nsim.le.5) then
     panres2@txFontHeightF = 0.024
  else
     panres2@txFontHeightF = 0.016
  end if
  if (SCALE_TIMESERIES.eq."True") then
     tt = ind(nyr.eq.nyr_max)
     panres2@gsnPanelScalePlotIndex = tt(0)
     delete(tt)
  end if
  if (nsim.le.12) then
     lp = (/nsim,1/)
  else
     lp = (/nrow,ncol/)   ;(/nsim/2+1,nsim/8+1/)  
  end if
  panres2@txString = "TAS Global Average (DJF)"
  gsn_panel2(wks_aa_djf,xy_djf,lp,panres2)  
  delete(wks_aa_djf)
  
  panres2@txString = "TAS Global Average (MAM)"
  gsn_panel2(wks_aa_mam,xy_mam,lp,panres2)  
  delete(wks_aa_mam)
  
  panres2@txString = "TAS Global Average (JJA)"
  gsn_panel2(wks_aa_jja,xy_jja,lp,panres2)  
  delete(wks_aa_jja)
  
  panres2@txString = "TAS Global Average (SON)"
  gsn_panel2(wks_aa_son,xy_son,lp,panres2)  
  delete(wks_aa_son)
  
  panres2@txString = "TAS Global Average (Annual)"
  gsn_panel2(wks_aa_ann,xy_ann,lp,panres2)  
  delete(wks_aa_ann)
  
  panres2@txString = "TAS Global Average (Monthly)"
  gsn_panel2(wks_aa_mon,xy_mon,lp,panres2)  
  delete(wks_aa_mon)
  
  panres2@txString = "TAS Running 8yr Trend (Monthly)"
  gsn_panel2(wks_rt_mon,xy_rt_mon(0,:),lp,panres2)  
  
  panres2@txString = "TAS Running 10yr Trend (Monthly)"
  gsn_panel2(wks_rt_mon,xy_rt_mon(1,:),lp,panres2)  
  
  panres2@txString = "TAS Running 12yr Trend (Monthly)"
  gsn_panel2(wks_rt_mon,xy_rt_mon(2,:),lp,panres2)  
  
  panres2@txString = "TAS Running 14yr Trend (Monthly)"
  gsn_panel2(wks_rt_mon,xy_rt_mon(3,:),lp,panres2)  
  
  panres2@txString = "TAS Running 16yr Trend (Monthly)"
  gsn_panel2(wks_rt_mon,xy_rt_mon(4,:),lp,panres2)  
  delete(wks_rt_mon)
  
  delete([/nrow,ncol,lp,map_djf,map_mam,map_jja,map_son,map_ann,map_mon,xy_djf,xy_mam,xy_jja,xy_son,xy_ann,xy_mon/])
  delete([/xy_rt_mon/])
  delete(panres2)
  if (isfilepresent2("obs_trefht")) then
     delete([/xy_obs_mam,xy_obs_jja,xy_obs_son,xy_obs_ann,xy_obs_mon/])
  end if
  OUTDIR = getenv("OUTDIR")
  if (wks_type.eq."png") then  
     do gg = 1,5
        if (isfilepresent2(OUTDIR+"tas.runtrend.mon.00000"+gg+".png")) then
           system("mv "+OUTDIR+"tas.runtrend.mon.00000"+gg+".png "+OUTDIR+"tas."+rt_nyr(gg-1)+"yr_runtrend.mon.png")
        end if
     end do
  else
     if (isfilepresent2(OUTDIR+"tas.runtrend.mon.ps")) then
        system("psplit "+OUTDIR+"tas.runtrend.mon.ps "+OUTDIR+"tas_rt")
        do gg = 1,5
           if (isfilepresent2(OUTDIR+"tas_rt000"+gg+".ps")) then
              system("mv "+OUTDIR+"tas_rt000"+gg+".ps "+OUTDIR+"tas."+rt_nyr(gg-1)+"yr_runtrend.mon.ps")
           end if
        end do
        system("rm "+OUTDIR+"tas.runtrend.mon.ps")
     end if
  end if
  delete(OUTDIR)
  print("Finished: tas.trends_timeseries.ncl")
end
