; Ver 1.1 2023-06-19
; Modified path due to run in lab computer
;

function shock_line, s2,s3,s4,p

  p=plot([s2,s2],[p.yrange[0],p.yrange[1]],linestyle='dot',color='crimson',/over)
  p=plot([s3,s3],[p.yrange[0],p.yrange[1]],linestyle='dot',color='crimson',/over)
  p=plot([s4,s4],[p.yrange[0],p.yrange[1]],linestyle='dot',color='crimson',/over)

end

pro smith_2001_ace_reproduce



  fm=file_search('D:\data\ace\mag\level_2_cdaweb\mfi_h3\*.cdf')  ; RTN magnetic field (?)
  fswe=file_search('D:\data\ace\swepam\level_2_cdaweb\swe_h0\2000\*.cdf')

  adata=read_cdf(fswe[0],VAR='*')
  aatt=read_cdf(fswe[0],VAR='*',/att)
  stop

  fswi=file_search('D:\data\ace\swics\level_2_cdaweb\swi_h6\2000\*.cdf')

  swevar=['V_RTN','NP','TPR','EPOCH']         ; NO DATA FOR LONGTIUDE, LATITUDE
  mvar=['MAGNITUDE','EPOCH','BGSM','BRTN']    ; NO DATA FOR LONGTIUDE, LATITUDE
  swivar=['EPOCH','NH','VH','VTHH']

  absB=[]
  Blat=[]
  Blon=[]
  BRTN=[]
  Btime=[]  ; Time resolution [1s]
  time=[]   ; Time resolution [1hr]
  Np=[]
  Tp=[]     ;
  Vr=[]     ;
  Bz=[]     ; in GSM coordinate
  SWIvr=[]
  SWInp=[]
  SWIVth=[]
  SWItime=[]


  for ii=0, n_elements(fm)-1 do begin

    ; MAGNITUDE MFI ---------------------------------------------------------------------
    adata=read_cdf(fm[ii],var=mvar)
    aatt=read_cdf(fm[ii],var=mvar,/att)

    ind=where(adata.MAGNITUDE EQ aatt.MAGNITUDE.FILLVAL, cnt)
    if cnt ne 0 then adata.MAGNITUDE[ind]=!values.f_nan
    ind=where(adata.BGSM EQ aatt.BGSM.FILLVAL, cnt)
    if cnt ne 0 then adata.BGSM[ind]=!values.f_nan
    ind=where(adata.BRTN EQ aatt.BRTN.fillval, cnt)
    if cnt ne 0 then adata.BRTN[ind]=!values.f_nan

    Bz=[Bz,adata.BGSM[*,2]]                     ; z component but southward?
    absB=[absB,adata.MAGNITUDE[*,1]]
    Btime=[Btime,cdf_epoch_tojuldays(adata.EPOCH[*,2])]
    BRTN=[BRTN,adata.BRTN]
    ;------------------------------------------------------------------------------------

    ; SWEPAM ----------------------------------------------------------------------------
    adata=read_cdf(fswe[ii],var=swevar)
    aatt=read_cdf(fswe[ii],var=swevar,/att)

    ind=where(adata.V_RTN EQ aatt.V_RTN.FILLVAL, cnt)
    if cnt ne 0 then adata.V_RTN[ind]=!values.f_nan
    ind=where(adata.NP EQ aatt.NP.FILLVAL, cnt)
    if cnt ne 0 then adata.NP[ind]=!values.f_nan
    ind=where(adata.TPR EQ aatt.TPR.FILLVAL, cnt)
    if cnt ne 0 then adata.TPR[ind]=!values.f_nan
    Np=[Np,adata.NP[*,0]]                               ; [time, three comp] but the same (check)
    Tp=[Tp,adata.TPR[*,0]]                              ; [time, three comp] but the same (check) [according to attribute, radial component of the proton temperature]
    Vr=[Vr,adata.V_RTN[*,0]]
    time=[time,cdf_epoch_tojuldays(adata.EPOCH[*,0])]   ; [time, three comp] but the same (check)
    ;------------------------------------------------------------------------------------

    ; SWICS ----------------------------------------------------------------------------
    adata=read_cdf(fswi[ii],var=swivar)
    aatt=read_cdf(fswi[ii],var=swivar,/att)

    ind=where(adata.NH EQ aatt.NH.FILLVAL, cnt)
    if cnt ne 0 then adata.NH[ind]=!values.f_nan
    ind=where(adata.VH EQ aatt.VH.FILLVAL, cnt)
    if cnt ne 0 then adata.VH[ind]=!values.f_nan
    ind=where(adata.VTHH EQ aatt.VTHH.FILLVAL, cnt)
    if cnt ne 0 then adata.VTHH[ind]=!values.f_nan

    SWIvr=[SWIvr, adata.VH[*,0]]
    SWInp=[SWInp, adata.NH[*,0]]
    SWIVth=[SWIVth, adata.VTHH[*,0]]
    SWItime=[SWItime, cdf_epoch_tojuldays(adata.EPOCH[*,0])]     ; time resolution = 12 min
  endfor

  ;------------------------------------------------------------------------------------
  ; Thermal temperature from thermal speed
  ;------------------------------------------------------------------------------------

  mp=1.67*1e-27     ; kg
  kb = 1.38*1e-23   ; K
  SWIVth*=1e3       ; km/s -> m/s (mks units)
  SWITp=mp*SWIVth^2/kb  ; 원래 2로 나눠야 하는 것 같은데 안 나눔....
  ;------------------------------------------------------------------------------------
  ;
  ;------------------------------------------------------------------------------------
  ; B lat, lon from RTN coordinate
  ; BRTN[*,0] = R component
  ; BRTN[*,1] = T component
  ; BRTN[*,2] = N component
  ;------------------------------------------------------------------------------------

  Blat=atan(BRTN[*,2]/sqrt(BRTN[*,0]^2+BRTN[*,1]^2))*180./3.141592

  Blon=make_array(n_elements(BRTN[*,0]),value=!values.f_nan)
  Bx=BRTN[*,0] ; R comp
  By=BRTN[*,1] ; T comp
  Bz=BRTN[*,2] ; N comp


  ind=where(Bx gt 0 and By gt 0, cnt)   ; 1사분면
  if cnt ne 0 then Blon[ind]=atan(By[ind]/Bx[ind])*180./3.141592
  ind=where(Bx gt 0 and By lt 0, cnt)   ; 4사분면
  if cnt ne 0 then Blon[ind]=(2.*!pi+atan(By[ind]/Bx[ind]))*180./3.141592
  ind=where(Bx lt 0 and By gt 0, cnt)   ; 2사분면
  if cnt ne 0 then Blon[ind]=acos(Bx[ind]/sqrt(Bx[ind]^2+By[ind]^2))*180./3.141592
  ind=where(Bx le 0 and By le 0, cnt)   ; 3사분면
  if cnt ne 0 then Blon[ind]=(!pi-asin(By[ind]/sqrt(Bx[ind]^2+By[ind]^2)))*180./3.141592


  ind=where(Bx gt 0. and By eq 0., cnt)     ; positive X-axis
  if cnt ne 0 then Blon[ind]=0.
  ind=where(Bx eq 0. and By gt 0., cnt)     ; positive Y-axis
  if cnt ne 0 then Blon[ind]=90.
  ind=where(Bx lt 0. and By eq 0., cnt)     ; negative X-axis
  if cnt ne 0 then Blon[ind]=180.
  ind=where(Bx eq 0. and By lt 0., cnt)     ; negative Y-axis
  if cnt ne 0 then Blon[ind]=270.
  ;------------------------------------------------------------------------------------

  ;  ;------------------------------------------------------------------------------------
  ;  ; Debugging
  ;  ;------------------------------------------------------------------------------------
  ;  p=plot(time, vr, xtickunit='day')
  ;  p=plot(switime, swivr, xtickunit='day',/over,color='red')
  ;
  ;
  ;
  ;  p=plot(time, Tp, xtickunit='day',/ylog)
  ;  p=plot(switime, swiTP, xtickunit='day',/over,color='red')
  ;  p=plot(time, Np, xtickunit='day')
  ;  p=plot(switime, swiNp, xtickunit='day',/over,color='red')

  ;------------------------------------------------------------------------------------
  ; smoothing
  ; 12 minute -> 1 hour
  ;------------------------------------------------------------------------------------

  new_time=[]
  new_Vr=[]
  new_Np=[]
  new_Tp=[]

  for tt=0, n_elements(switime)-3 do begin
    new_time=[new_time, mean(switime[tt:tt+2],/NaN)]
    new_Vr=[new_Vr, mean(swiVr[tt:tt+2],/NaN)]
    new_Np=[new_Np, mean(swiNp[tt:tt+2],/NaN)]
    new_Tp=[new_Tp, mean(swiTp[tt:tt+2],/NaN)]
  endfor

  ;------------------------------------------------------------------------------------
  ; To replace SWEPAM with SWICS data (VER 2)
  ;------------------------------------------------------------------------------------
  st=julday(01,196,2000,11,00,00)
  ed=julday(01,198,2000,01,00,00)

  swi_ind=where(new_time ge st and new_time le ed, cnt)

  SWItime = new_time[swi_ind]
  SWIVr =new_Vr[swi_ind]
  SWINp=new_Np[swi_ind]
  SWITp=new_Tp[swi_ind]
  ;------------------------------------------------------------------------------------


  ;------------------------------------------------------------------------------------
  ; To replace SWEPAM with SWICS data (VER 1)
  ;------------------------------------------------------------------------------------
  ;  st=julday(01,196,2000,11,00,00)
  ;  ed=julday(01,198,2000,01,00,00)
  ;
  ;  swi_ind=where(switime ge st and switime le ed, cnt)
  ;
  ;  SWItime = switime[swi_ind]
  ;  SWIVr =SWIvr[swi_ind]
  ;  SWINp=SWINp[swi_ind]
  ;  SWITp=SWITp[swi_ind]
  ;------------------------------------------------------------------------------------


  ;------------------------------------------------------------------------------------
  ; Plasma beta & Alfven speed
  ; to adjust time resolution
  ;------------------------------------------------------------------------------------
  ; FOR SWEPAM
  time_array=value_locate(Btime,time)
  kb = 1.38*1e-23
  Npm=Np*1e6          ; cm^-3 -> m^-3
  absBT=absB[time_array]*1e-9
  betaa =  8.*!pi*1e-7*Npm*kb*Tp/absBT/absBT

  muu=4.*!pi*1e-7
  mp=1.67*1e-27     ;kg
  alfven_v=sqrt(absBt*absBt/muu/mp/Npm)/1e3     ;m/s->km/s

  ; FOR SWICS
  time_array=value_locate(Btime,SWItime)
  Npm=SWINp*1e6          ; cm^-3 -> m^-3
  absBT=absB[time_array]*1e-9
  SWIbetaa =  8.*!pi*1e-7*Npm*kb*SWITp/absBT/absBT
  SWIalfven_v=sqrt(absBt*absBt/muu/mp/Npm)/1e3     ;m/s->km/s
  ;------------------------------------------------------------------------------------


  w2=window(dimension=[600,800])
  ;w2=window(dimension=[1200,1800])
  pos=make_plot_pos(0.12,0.05,0.95,0.95,9,gap=0.005,/col)
  ld=label_date(date_format='%N/%D')

  ex={xtickformat:'label_date',xminor:1,xstyle:1,current:w2,font_size:9,symbol:'dot',sym_size:0.5,xrange:[julday(07,12,2000,0,0,0),julday(07,18,2000,0,0,0)],xtickvalues:[julday(07,12,2000,0,0,0),julday(07,13,2000,0,0,0),julday(07,14,2000,0,0,0),julday(07,15,2000,0,0,0),julday(07,16,2000,0,0,0),julday(07,17,2000,0,0,0),julday(07,18,2000,0,0,0)]}
  tt=text(0.5,0.97,'ACE observations',alignment=0.5,/current)

  ; SHOCK
  s2=julday(01,195,2000,09,19,13)
  s3=julday(01,196,2000,14,59,49)
  s4=julday(01,197,2000,14,16,24)


  p=scatterplot(Btime,absB,_extra=ex,pos=pos[*,-1],ytitle='B [nT]',yminor=1) ; absB
  p.xshowtext=0
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(Btime,Blat, yrange=[-90.,90],_extra=ex,pos=pos[*,-2],ytitle='B$_l_a_t, \delta$',ytickinterval=45,yminor=0)
  p.xshowtext=0
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(Btime, Blon, yrange=[0,360],_extra=ex,pos=pos[*,-3],ytitle='B$_l_o_n, \lambda$',ytickinterval=90,yminor=0); Blon
  p.xshowtext=0
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(Btime,Bz,_extra=ex,pos=pos[*,-4],ytitle='$B_s$(GSM)',yrange=[-60,30],ytickinterval=20,yminor=0)     ;B_GSM_southward
  p.xshowtext=0
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(time,Vr,_extra=ex,pos=pos[*,-5],ytitle='$V_R$',yrange=[300,1100],ytickvalues=[400,600,800,1000],yminor=0)
  p.xshowtext=0
  p=plot(SWItime, SWIVr, /over, color='red')
  tmp=shock_line(s2,s3,s4,p)


  p=scatterplot(time,Np,_extra=ex,pos=pos[*,-6],ytitle='$N_P [cm^-^3]$',/ylog,ytickunit='scientific',ytickvalue=[0.1,1,10])
  p.xshowtext=0
  p=plot(SWItime, SWINp, /over, color='red')
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(time,Tp,_extra=ex,pos=pos[*,-7],ytitle='$T_P [K]$',/ylog,yminor=0,ytickunit='scientific',ytickvalue=[1e4,1e5])
  p.xshowtext=0
  p=plot(SWItime, SWITp, /over, color='red')
  tmp=shock_line(s2,s3,s4,p)

  p=scatterplot(time,betaa,_extra=ex,pos=pos[*,-8],ytitle='$\beta$',/ylog,yrange=[5*1e-4,1e1],ytickunit='scientific',yminor=0,ytickvalue=[0.001,0.01,0.1,1])
  p.xshowtext=0
  p=plot(SWItime, SWIbetaa, /over, color='red')
  tmp=shock_line(s2,s3,s4,p)


  p=scatterplot(time,Alfven_v,_extra=ex,pos=pos[*,-9],ytitle='$V_A [km/s]$',xtitle='Day',yrange=[0,1000],ytickvalues=[0,400,800])
  p=plot(SWItime, SWIAlfven_v, /over, color='red')
  tmp=shock_line(s2,s3,s4,p)

  stop
  stop
  ;p.save,'C:\Users\neagg\OneDrive\바탕 화면\task\lab\temp\20230614\smith_2001_reproduce1_ACE_scatterplot.png'
  ;w2.save,'E:\temp\20230619\smith_2001_reproduce1_ACE_scatterplot.png'
  w2.save,'D:\temp\20230620\smith_2001_Figure1_ACE_scatterplot_SWICS_12min.png'
  w2.save,'D:\temp\20230620\smith_2001_Figure1_ACE_scatterplot_SWICS_30minsmoothing.png'
  stop
  stop


  stop
  stop



end

;------------------------------------------------------------------------------------
; DEBUGGING
;------------------------------------------------------------------------------------
;p=plot(SWIvr[*,0])
;p=plot(SWIvr[*,1],/over,color='red')
;p=plot(SWIvr[*,2],/over,color='blue')
;
;p=plot(SWINp[*,0])
;p=plot(SWINp[*,1],/over,color='red')
;p=plot(SWINp[*,2],/over,color='blue')
;
;p=plot(SWITp[*,0])
;p=plot(SWITp[*,1],/over,color='red')
;p=plot(SWITp[*,2],/over,color='blue')
;
;p=plot(SWITime[*,0])
;p=plot(SWITime[*,1],/over,color='red')
;p=plot(SWITime[*,2],/over,color='blue')
;------------------------------------------------------------------------------------


;  ;------------------------------------------------------------------------------------
;  ; Debugging
;  ;------------------------------------------------------------------------------------
;  Bx=1. & By=sqrt(3)
;  ind=where(Bx ge 0 and By ge 0, cnt)   ; 1사분면
;  print, cnt
;  print, atan(By/Bx)*180./3.141592
;
;  Bx=1. & By=-sqrt(3)
;  ind=where(Bx ge 0 and By le 0, cnt)   ; 4사분면
;  print, cnt
;  print, (2.*!pi+atan(By/Bx))*180./3.141592
;
;  Bx=-1. & By=sqrt(3.)
;  ind=where(Bx le 0 and By ge 0, cnt)   ; 2사분면
;  print, cnt
;  print, acos(Bx/sqrt(Bx^2+By^2))*180./3.141592
;
;  Bx=-1. & By=-sqrt(3.)
;  ind=where(Bx le 0 and By le 0, cnt)   ; 3사분면
;  print, (!pi-asin(By/sqrt(Bx^2+By^2)))*180./3.141592
;  ;------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------
; Debugging
;------------------------------------------------------------------------------------
;  for ii=0, n_elements(ind2)-1 do begin
;    ind_what=where(ind1 eq ind2[ii], cnt)
;    if cnt ne 0 then stop
;  endfor
;
;  for ii=0, n_elements(ind4)-1 do begin
;    ind_what=where(ind3 eq ind4[ii], cnt)
;    if cnt ne 0 then stop
;  endfor
;------------------------------------------------------------------------------------
;
;fswe=file_search('C:\Users\neagg\data\ace\swepam\swe_h2\2000\*.cdf')
;
;pos_RTN=adata.RTN_ATT
;time=conv_time(adata.epoch, /from_epoch, /to_julday)
;
;p=plot(time[*,0])
;p=plot(time[*,1],/over,color='blue')
;p=plot(time[*,2],/over,color='red')
;
;adata=read_cdf(fm[0],VAR=mvar)
;aatt=read_cdf(fm[0],VAR=mvar,/att)
;
;ind=where(adata.MAGNITUDE LE)
;
;
;adata=read_cdf(fswe[0],VAR=swevar)
;aatt=read_cdf(fswe[0],VAR=swevar,/att)


;  fm=file_search('C:\Users\neagg\data\ace\mag\level_2_cdaweb\mif_h2\2000\*.cdf')
;  fm=file_search('C:\Users\neagg\data\ace\mag\level_2_cdaweb\mfi_h0\2000\*.cdf')


;  adata=read_cdf(fswi[0],VAR=swivar)
;  aatt=read_cdf(fswi[0],VAR=swivar,/att)




