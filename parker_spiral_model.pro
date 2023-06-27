;Ver. 1. 2023-06-21

pro parker_spiral_model

  ;------------------------------------------------------------------------
  ; constant
  ;------------------------------------------------------------------------
  r0  = 0.0001                        ; AU     (solar surface)
  w   = 2.7*1e-6                      ; rad/s  (rotation speed)
  ;bd=[-1.5,1.5]                           ; boundary (display boundary)
   bd=[-5,5]                           ; boundary (display boundary)


  ;------------------------------------------------------------------------
  ; setting
  ;------------------------------------------------------------------------
  phi0_deg = [30.,120.,210.,300.]
  kmtoau   = 1.496*1e8               ; 1 AU = 1.496*1e8 km
  v0       = 700./kmtoau             ; km/s -> AU/s         (solar wind speed)
  phi0 = phi0_deg*3.141592/180.      ; convert deg to rad
  r = linspace(0.,5.,1000)        ; AU
  disR = r-r0                        ; AU                   (r-r0)
  cr = ['crimson','green','orange','sky blue']

  ;------------------------------------------------------------------------
  ; frame
  ;------------------------------------------------------------------------
  ww=window(dimension=[400,400])
  pos=make_plot_pos(0.15,0.18,0.80,0.88,1,/col)
  exlg = {shadow:0,font_size:11,auto_text_color:1,sample_width:0}
  p=plot(bd,[0,0],color='gray',xrange=bd,yrange=bd,title='SW speed = '+strcompress(string(v0*kmtoau, format='(F4.0)'),/rem)+' km/s',pos=pos,xtitle='X [AU]',ytitle='Y [AU]',font_name='Calibri',font_size=11,/current)
  p=plot([0,0],bd,color='gray',/over)

  xx=1.*cos(linspace(0.,2.*!pi,30))
  yy=1.*sin(linspace(0.,2.*!pi,30))

  p=plot(xx,yy,linestyle='--',thick=1.5,color='gray',/over)
  lgd=[]

  for phi_ind=0, n_elements(phi0_deg)-1 do begin
    ;------------------------------------------------------------------------
    ; Equation
    ;------------------------------------------------------------------------
    phi     = phi0[phi_ind] - w*disR/v0         ; Radian
    phi_deg = phi*3.141592/180.                 ; Deg
    ;------------------------------------------------------------------------

    xx = r*cos(phi)
    yy = r*sin(phi)

    p=plot(xx,yy,/over,thick=2, color=cr[phi_ind],name=strcompress(string(phi0_deg[phi_ind],format='(I3)'),/rem)+'$\deg$')
    lgd=[lgd,p]
  endfor
  lg=legend(target=lgd,_extra=exlg,position=[0.85,0.5],vertical_alignment=0.5,horizontal_alignment=0,font_name='calibri')
  stop
  p.save,'D:\temp\20230621\parker_spiral_model_'+strcompress(string(v0*kmtoau, format='(I3)'),/rem)+'_zoom_in.png'
  p.save,'D:\temp\20230621\parker_spiral_model_'+strcompress(string(v0*kmtoau, format='(I3)'),/rem)+'.png'

  stop

end