external CLASSIFY "./algorithm_class_precip_morphology.so"

;define path to WRF output
path = "/glade/scratch/rberrios/timevaryingshear_wrf/run_fromctl_3ms_todeep_15ms_at"+sprinti("%0.3i",fhr_shr)+"_ndg12hfull/"

;get WRF list of files
files = systemfunc("ls -1 "+path+"wrfout_d03*00:00")

;define an array with the forecast hours (assumed to be hourly output)
;this is in case we want to choose a particular time interval (eg., 6-hly)
fhrvec = ispan(0,dimsizes(files)-1,1)
nhrs = dimsizes(fhrvec)

;define radii for area-averaging the number of convective, stratiform, and anvil points
R = (/100,150,200/)
nr = dimsizes(R)

;define variable that will hold the number of convective, stratiform, and anvil
 points per each forecast time
pmorph = new((/nhrs,nr,4/),integer) ;4 classes
;define variable that will hold the number of deep, moderate, and shallow convective points per each forecast time
conv_type = new((/nhrs,nr,3/),integer) ;3 types of convection

;get center info
reflev=900
file_with_ctr_info = path+"/hourly_z_object_centroid_ind_"+reflev+"hPa_d03"
data = asciiread(file_with_ctr_info,-1,"string")
delim = " "
yposmid = str_get_field(data,1,delim)
xposmid = str_get_field(data,2,delim)
yposmid@_FillValue = -2147483647
xposmid@_FillValue = -2147483647
xposmid := stringtointeger(xposmid(1::));remove header
yposmid := stringtointeger(yposmid(1::))
;apply a 6-h running mean
yposmid := runave(yposmid,6,6)
xposmid := runave(xposmid,6,6)

;we're ready to loop over the WRF output files
do hh = 0, nhrs-1
  fhr = fhrvec(hh)
  print("working on "+fhr+"...")
  ;read reflectivity and w, and interpolate to 0.9 km and 3 km
  ncid = addfile(files(fhr),"r")
  refall = wrf_user_getvar(ncid,"dbz",0)
  z = wrf_user_getvar(ncid,"z",0)
  var2plt = wrf_user_intrp3d(refall,z,"h",(/400.0,3000./),0.,False)
  www = wrf_user_intrp3d(wrf_user_getvar(ncid,"wa",0),z,"h",(/400.0,3000./),0.,False)
  var2plt@_FillValue = -9999.
  ;read in grid spacing info
  DX = ncid@DX*1.e-3
  DY = ncid@DY*1.e-3

  ;call in external FORTRAN function to do the classification
  dims = dimsizes(var2plt)
  varout = new(dims(1::),float)
  CLASSIFY :: class_precip_morphology(dims(2),dims(1),dims(0),var2plt,www,varout,DY,DX)

  ;---- added on 12/7/2020: for convective points, find level of 20 dBZ echo top ---------;
  varout_conf = conform(z,varout,(/1,2/))
  ; First pass: find convective points
  refall_conv = where(varout_conf .eq. 1., refall, refall@_FillValue)
  ; Second pass: find the 20dBZ contour
  refall_conv = where(refall_conv .ge. 17.5 .and. refall_conv .le. 22.5,refall_conv,refall@_FillValue)
  ; Last pass: find the height of the 20-dBZ contour? 
  z_conv_20dbZ = dim_max_n(where(.not. ismissing(refall_conv), z, z@_FillValue),0)
  varout_conv = new(dims(1::),integer)
  ;shallow conv
  varout_conv = where( .not. ismissing(z_conv_20dbZ) .and. z_conv_20dbZ .le. 6000.0,1,varout_conv)
  ;moderate conv
  varout_conv = where( .not. ismissing(z_conv_20dbZ) .and. z_conv_20dbZ .ge. 6000.0 .and. z_conv_20dbZ .le. 10000.0,2,varout_conv)  
  ;deep conv
  varout_conv = where( .not. ismissing(z_conv_20dbZ) .and. z_conv_20dbZ .ge. 10000.0,3,varout_conv) 
  delete(varout_conf)
  delete(refall_conv) 

  ;define variables that we'll use to zero out domain surrounding outside the desired radii
  ictr = yposmid(fhr)
  jctr = xposmid(fhr)

  ydum = (ispan(0,dims(1)-1,1)-ictr)*DY
  xdum = (ispan(0,dims(2)-1,1)-jctr)*DX
  y2D = conform_dims(dims(1::),ydum,0)
  x2D = conform_dims(dims(1::),xdum,1)
  dist = sqrt(x2D*x2D+y2D*y2D)

  ;loop over desired radii
  do rr = 0, nr-1
    var2save := where(dist .gt. R(rr),varout@_FillValue,varout)
    ;Now get only points within area of interest
    var1D := ndtooned(var2save)
    indAI := ind(.not.ismissing(var1D))
    if (.not. any(ismissing(indAI)))
      var1DAI := var1D(indAI)
      pmorph(fhr,rr,:) = (/dimsizes(ind(var1DAI.eq.1.)),dimsizes(ind(var1DAI.eq.2.)),dimsizes(ind(var1DAI.eq.3.)),dimsizes(ind(var1DAI.eq.0.))/)
    end if

    ;repeat for shallow, moderate, and deep convection
    var2save := where(dist .gt. R(rr),varout_conv@_FillValue,varout_conv)
    ;Now get only points within area of interest
    var1D := ndtooned(var2save)
    indAI := ind(.not.ismissing(var1D))
    if (.not. any(ismissing(indAI)))
      var1DAI := var1D(indAI)
      conv_type(fhr,rr,:) = (/dimsizes(ind(var1DAI.eq.1.)),dimsizes(ind(var1DAI.eq.2.)),dimsizes(ind(var1DAI.eq.3.))/)
    end if
  end do ;--- end of loop over radii
end do ;--- end of loop over output files

  ; plot reflectivity and classification to make sure everything works
  wks = gsn_open_wks("png","example_pmorph_"+fhr)
  gsn_define_colormap(wks,"radar_DuPage")

  res = True
  res@gsnFrame = False
  res@gsnDraw = False
  res@cnFillOn = True
  res@cnLinesOn = False
  res@tiXAxisString = "distance (km)"
  res@tiYAxisString = "distance (km)"
  res@pmLabelBarOrthogonalPosF = 0.1
  res@tiXAxisOffsetYF = 0.1
  res@pmLabelBarHeightF = 0.05
  ;res@lbTitleOn = True
  res@lbTitlePosition = "Bottom"
  res@lbTitleFontHeightF = 0.012
  res@gsnCenterString = ""
  res@gsnRightString = ""
  res@gsnLeftString = ""
  res@cnLevelSelectionMode = "ExplicitLevels"
  res@cnLevels = ispan(5,75,5)
  res@cnFillPalette = "radar_DuPage"
  ;res@gsnLeftString  = "t = "+fhr+" h"
  ;res@gsnRightString = "N = "+(cc)
  res@vpWidthF  = 0.38
  res@vpHeightF = 0.38
  res@vpXF = 0.13
  res@vpYF = 0.54
  res@lbTitleString = "0.9-km reflectivity (dBZ)"
  ;res@trXMinF = 180-25
  ;res@trXMaxF = 180+25
  ;res@trYMinF = 180-25
  ;res@trYMaxF = 180+25
  res@trYMinF = -250
  res@trXMinF = -250
  res@trYMaxF = 250
  res@trXMaxF = 250
  res@cnFillMode = "RasterFill"

  var2cnt = var2plt(0,:,:)
  var2cnt!1 = "x"
  var2cnt!0 = "y"
  var2cnt&x = (ispan(1,dims(2),1)-dims(2)/2.)*2.
  var2cnt&y = (ispan(1,dims(1),1)-dims(1)/2.)*2.
  copy_VarCoords(var2cnt,varout)
  plot1 = gsn_csm_contour(wks,var2cnt,res)


  res2 = res
  res2@vpXF = 0.58
  res2@tiYAxisOn = False
  res2@cnLevels := (/0.5,1.5,2.5/)
  res2@cnFillPalette = "SVG_Lindaa07"
  res2@cnFillColors = (/-1,219,80,191/)
  res2@cnExplicitLabelBarLabelsOn = True
  res2@lbLabelStrings = (/"no rain","convective","stratiform","other"/)
  res2@lbLabelAlignment = "BoxCenters"
  res2@lbTitleString = "precipitation"
  plot2 = gsn_csm_contour(wks,varout,res2)

  draw(plot1)
  draw(plot2)

  frame(wks)

end do ;--- end of loop over output files

; prepare to save a netcdf for analysis and plotting later on
pmorph!0 = "time"
pmorph!1 = "R"
pmorph!2 = "class"
pmorph&time = fhrvec
pmorph&R = R
pmorph&class = (/1,2,3,0/)
pmorph@class_description = (/"1:convective","2:stratiform","3:other","0:no_rain"/)

conv_type!0 = "time"
conv_type!1 = "R"
conv_type!2 = "height"
conv_type&height = (/1,2,3/)
conv_type@height_description = (/"1: shallow","2:moderate","3:deep"/)

;save netCDF file
ncout = addfile(path+"precip_morph_radii_"+reflev+"hPactr.nc","c")
ncout->pmorph = pmorph
ncout->conv_type=conv_type
delete(ncout)
