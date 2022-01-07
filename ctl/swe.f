      program SILEPE                                                             
c                                                                               
c                                                                             
c  This  program  performs the single  level  primitive          
c  equation forecast. It can also be used for iterative          
c  initialization.  The imitive equations use the zonal          
c  and  meridional  components  of  the wind  field and          
c  the height field  as the predicted variables.          
c  The advective  terms  are  computed  using  a  semi-          
c  lagrangian approach.time integration is accomplished
c  through the  Matsuno time  integration scheme.                                                       
c  The  model  can  be used with or without the terrain          
c  field  and  its associated  gradients.  these  field          
c  are assumed to be  given. the model also allows  for          
c  closed or open boundaries. the domain invariants are          
c  computed during each output time .                             
c                                                                             
c                                                                             
c  originator     : Florida State University,last revision
c                   L. Bounoua, 1994                                     
c                                                                             
c  input          : (1) logical  unit 21 contains the terrain  field           
c                   its  associated zonal and meridional  gradients.          
c                   each array of data  is stored in one  record and          
c                   are read in the  zm,  dhdx,  dhdy  order .          
c               (2) logical  unit 22  contains the initial zonal (u)          
c                   and meridional (v) components  of the wind field          
c                   and height field (z).  each array  is stored  in          
c                   one  record in the  u,v,z order.          
c               (3) logical unit 23 contains information required by          
c                   the program during execution.  see comment cards          
c                   below for details.                                        
c                                                                             
c  output         : (1) logical unit 28 contains the output  fields.          
c                   the first 3 records contains the initial  u,v,z.          
c                   subsequent  sets   of  3   records  contain  the          
c                   forecast u,v,z at 12 hourly intervals ( or other          
c                   intervals that has been set by the user.)                 
c               (2) logical unit 6 is used for listing.                      
c                                                                             
c       
c                                                                               
c  Variables input through unit 23.        
c                                                                               
c  variable name   description                    read in 
c  -------------   -----------                    -------            
c  ncase           description of forecast         const              
c  time            starting time of forecast       const              
c  ijk  (=0)       iterative initialization        const              
c       (=1)       one level pe forecast                              
c  iterr   (=0)    no terrain effects included     const              
c          (=1)    terrain included                                   
c  ioutpt  (=0)    no gain and shade done          const              
c          (=1)    gain and shade done                                
c  nsbd    (=0)    v=0. at n-s boundary            const              
c          (=1)    smoothing done at n-s boundary  const              
c  alfa            smoothing coeff (0. - 1.)       const              
c  l,m,dphi        zonal dimension, meridional                        
c                  dimension and grid size in deg  const              
c  slat,dt         southern most lat, time step    const              
c  nohist          history tape parameter          const              
c                  (=999)  no history tapes made                              
c
c  iterative initialization                                                         
c  -------------------------                                                          
c  niter           number of iterations            const              
c  nout            number of iter. before write    const              
c                  (=1 for one time step)                             
c  itcnt           no. of history rec. to be read  const              
c
c  forecast                                                           
c  --------                                                           
c  tfcst           length of forecast in hours     const              
c  toint           output intervals                const              
c   tcnt           no. of sets of records to be    const              
c                  read                                               
c                                                                               
      parameter (lg = 114,mg = 17)                                                    
c                                                                               
      common /a01/  l,      m,       l1,      m1,       m2,      g,             
     &              slat,   dphi                                                
      common /a02/  jin,    jout,    itape,   jb,       jc,      jd             
      common /a03/  ijk,    itcnt,   itime,   nstep,   iterr                    
      common /a04/  time,   dt,      nohist,  ioutpt,   nsbd,    alfa           
      common /a05/  tfcst,  toint,   niter,   nout,     title(3)                
      common /a06/  denom(6,mg), dx(mg),       cor(mg),    dy                   
      common /a07/  dx2iv(mg),   dxsqiv(mg),   dy2iv,      dysqiv               
      common /a08/  zm(lg,mg),   dhdx(lg,mg),  dhdy(lg,mg)                      
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg)                         
      common /a10/  up(lg,mg),   vp(lg,mg),    zp(lg,mg)                        
      common /a11/  uo(lg,mg),   vo(lg,mg),    zo(lg,mg)                        
      common /a12/  a(lg,mg),    b(lg,mg),     c(lg,mg)                         
      integer irec,debug,restart,iday,ihour,imin
      character*100 ofile
      irec    = 1
      restart = 1
      ofile   = 'ctl_00:00:00.dat'
c
c  Open input and output files .
c
      iterr   = 0
      if(iterr.ne.0) then                                                    
      open (21,file='terrain.out    ',status='old')
      endif                                  
      open (22,file='bgd.dat',status='old')                              
      open (28,file='ctl.dat',form='unformatted',
     & access='direct',recl=lg*mg*4)
c                                                                               
c  initialize variables in the common block.                                                                           
c                                                                               
      call ZEROS                                                                 
c                                                                               
c  Call subroutine CONST to define constants               
c  parameters and constants required by the program              
c                                                                               
      call CONST2(debug)                                                                 
c                                                                               
      if (debug.ge.1) then
      write(6,1000) time                                                          
 1000 format(5x,'starting time of forecast is ',f5.0)                           
 1001 format(30x,'iterative initialization')                                  
 1002 format(10x,'number of iterations =',i3,2x,'no of iterations ',          
     &'between outputs =',i2,2x,'no of history records to be read =',         
     & i3)                                                                     
 1003 format(30x,'  pe. forecast ')                                           
 1004 format(10x,'length of forecast in hrs =',f10.3,3x,                       
     &' output interval in hrs =',f6.2,3x,
     &10x,' no of history records,to be read =',i5)                                                       
 1005 format(10x,'l=',i5,2x,'m=',i5,2x,'slat=',f5.1,2x,                         
     &'dphi=',f6.3,2x,'dt=',f5.1)                                              
c                                                                               
      if (ijk.eq.0) then                                                          
      write(6,1001)                                                            
      write(6,1002) niter, nout, itcnt                                          
      else                                                                      
  45  continue                                                                
      write(6,1003)                                                             
      write(6,1004) tfcst, toint, itcnt                                         
      endif                                                                     
      write(6,1005) l, m, slat, dphi, dt
      endif
c                                                                               
c  Call subroutine INDATA to  read in the               
c  terrain and initial fields.             
c                                                                               
      call INDATA
      open(29,file=ofile,status='unknown')
      write(29,'(6F13.6)')((u(i,j),i=1,lg),j=1,mg)
      write(29,'(6F13.6)')((v(i,j),i=1,lg),j=1,mg)
      write(29,'(6F13.6)')((z(i,j),i=1,lg),j=1,mg)
      close(29)
      write(28,rec=irec)((u(i,j),i=1,lg),j=1,mg)
      irec    = irec + 1
      write(28,rec=irec)((v(i,j),i=1,lg),j=1,mg)
      irec    = irec + 1
      write(28,rec=irec)((z(i,j),i=1,lg),j=1,mg)
      irec    = irec + 1
c
c  Call subroutine INVART to compute the
c  closed domain invariants.
c                                                                
      call INVART(debug)                                                                
c                                                                               
      do 11100 j= 1, m,m1                                                          
      do 11100 i= 1, l                                                             
         uo(i,j)= u(i,j)                                                        
         vo(i,j)= v(i,j)                                                        
11100    zo(i,j)= z(i,j)                                                        
      jin       = ifix((toint*3600.)/dt)*itcnt                                          
      if (itcnt.eq.1) jin = 0                                                      
      if (ijk.eq.0  ) then                                                         
      jint      = nout*nstep                                                       
      if (jint.eq.0 ) jint = 1                                                  
      jout      = jin + jint                                                       
      stop      = niter/nout                                                      
      else                                                                      
      jint      = toint*3600./abs(dt)                                              
      if (jint.lt.1) jint = 1                                                    
      jout      = jin+jint                                                           
      jstop     = tfcst/toint+0.5                                                   
      endif                                                                     
c                                                                               
c  The 'do loop 11102' controls  the outputting  of the                
c  forecast    products.   during    each    looping,                
c  subroutine  'FCST'  is  first  called  to  perform                
c  the  time  integration   until  the   output  time                
c  is reached. the forecast products are  then output                
c  to  the output file if  required.  Next subroutine                
c  'INVART' is called  to compute  the closed  domain                
c  invariants of the model.                                          
c                                                                               
      do 11102 kk= 1, jstop                                                         
        print*,kk,time,jin,jout,jint,toint
        call FCST(debug)                                                                
        write(28,rec=irec)((u(i,j),i=1,lg),j=1,mg)
        irec    = irec + 1
        write(28,rec=irec)((v(i,j),i=1,lg),j=1,mg)
        irec    = irec + 1
        write(28,rec=irec)((z(i,j),i=1,lg),j=1,mg)
        irec    = irec + 1
        call INVART(debug)                                                              
        if (debug.ge.1) then
         write(6,1006) kk                                                        
 1006    format(2x,'after invart',5x,'output step number = ',i5)
        endif 
        iday      = ifix (time/86400.)
        ihour     = ifix ((time-iday*86400.)/3600.)
        imin      = ifix ((time-iday*86400.-ihour*3600.)/60.)
        ofile   = 'ctl_00:00:00.dat'
        if (iday.lt.10) then
         write(ofile(6:6),'(1I1)')iday
        elseif (iday.lt.100) then
         write(ofile(5:6),'(1I2)')iday
        else
         print*,'truth.exe: day string length is too long...stop'
         stop
        endif
        if (ihour.lt.10) then
         write(ofile(9:9),'(1I1)')ihour
        elseif (iday.lt.100) then
         write(ofile(8:9),'(1I2)')ihour
        else
         print*,'truth.exe: hour string length is too long...stop'
         stop
        endif
        if (imin.lt.10) then
         write(ofile(12:12),'(1I1)')imin
        elseif (iday.lt.100) then
         write(ofile(11:12),'(1I2)')imin
        else
         print*,'truth.exe: minute string length is too long...stop'
         stop
        endif 
        if (debug.ge.1) print*,'truth.exe: open truth output is: ',ofile(1:20)
        open(29,file=ofile,status='unknown')
        write(29,'(6F13.6)')((u(i,j),i=1,lg),j=1,mg)
        write(29,'(6F13.6)')((v(i,j),i=1,lg),j=1,mg)
        write(29,'(6F13.6)')((z(i,j),i=1,lg),j=1,mg)
        close(29)
        jin     = jout                                                              
        jout    = jin+jint                                                         
11102 continue                                                                  
      end
c                                                                       
      subroutine CONST2(debug)                                                           
c                                                                               
c  This routine computes or defines the constants that           
c  are needed in the program.  the values are  read in           
c  from logical unit 23 using free format read.                  
c                                                                               
      parameter(lg = 114,mg = 17)                                                    
c                                                                               
      common /a01/  l,      m,       l1,      m1,       m2,      g,             
     1              slat,   dphi                                                
      common /a02/  jin,    jout,    itape,   jb,       jc,      jd             
      common /a03/  ijk,    itcnt,   itime,   nstep,   iterr                    
      common /a04/  time,   dt,      nohist,  ioutpt,   nsbd,    alfa           
      common /a05/  tfcst,  toint,   niter,   nout,     title(3)                
      common /a06/  denom(6,mg), dx(mg),       cor(mg),     dy                  
      common /a07/  dx2iv(mg),   dxsqiv(mg),   dy2iv,       dysqiv 
      real      dy
      integer   debug
c                                                                               
c  reading namelist and set some internal values that were
c  previously input from tape23.dat. 
c
      call input_namelist(debug,l,m,tfcst,slat,dt,dy,toint)
      dphi      = dy/111e3
      ijk       = 1
      nohist    = 1
      times     = 0
      iterr     = 0
      ioutpt    = 0
      niter     = 5
      nout      = 5
      itcnt     = 1
      nsbd      = 1
      if (debug.eq.1) then
c       write (* ,*)l,m,ijk                                                        
c       write (* ,*)times,tfcst,toint,dt                                           
c       write (* ,*)iterr,ioutpt,niter,nout,itcnt,nsbd,nohist                      
c       write (* ,*)slat,dphi                                                      
c       write (* ,*)
       print*,'swe.exe: domain size l          =',l
       print*,'swe.exe: domain size m          =',m
       print*,'swe.exe: forecast time tfcst    =',tfcst
       print*,'swe.exe: time step dt           =',dt
       print*,'swe.exe: restart toint          =',toint
       print*,'swe.exe: domain grid dx         =',dy,dphi
       print*,'swe.exe: starting lat slat      =',slat
       read*
      endif 
c                                                                               
      m1        = m-1                                                                  
      l1        = l-1                                                                  
      m2        = m-2                                                                  
      time      = times*3600.                                                        
      nstep     = 6                                                                 
      alfa      = 0.95                                                               
      g         = 9.81                                                                  
      pi        = 4.*atan(1.)                                                          
      rad       = pi/180.                                                             
      dy        = 111.1*dphi*1.e03                                                    
      dy2iv     = 1./(2.*dy)                                                          
      dysqiv    = 1./(dy*dy)                                                         
      do 11110 j= 1, m                                                               
         x      = (slat+(j-1)*dphi)*rad                                               
         dx(j)  = dy*cos(x)                                                       
       dx2iv(j) = 1./(2.*dx(j))                                                  
      dxsqiv(j) = 1./(dx(j)*dx(j))                                            
        cor(j)  = 2.*7.292e-05*sin(x)                                            
11110 continue                                                                  
c                                                                               
c  denom is an array of interpolation constants 
c  used in the semi-lagrangian advective scheme           
c                                                                               
      dyy       = dy                                                                  
      do 11111 j= 2, m1                                                              
      denom(1,j)=  1./((2.*dx(j-1)*dx(j-1))*(2.*dyy*dyy))                    
      denom(2,j)= -1./((   dx(j-1)*dx(j-1))*(2.*dyy*dyy))                    
      denom(3,j)= -1./((2.*dx(j  )*dx(j  ))*(   dyy*dyy))                    
      denom(4,j)=  1./((   dx(j  )*dx(j  ))*(   dyy*dyy))                    
      denom(5,j)=  1./((2.*dx(j+1)*dx(j+1))*(2.*dyy*dyy))                    
      denom(6,j)= -1./((   dx(j+1)*dx(j+1))*(2.*dyy*dyy))                    
11111  continue                                                                  
      return                                                                    
      end                                                                       
c                                                                                
      SUBROUTINE input_namelist(debug,nx,ny,tfcst,
     &           slat,dt,dy,restart)
      INCLUDE "../registry/swe.inc"
      RETURN
      END

      subroutine FCST(debug)                                                           
c                                                                               
c  This  subroutine  does the time integration part of           
c  the model.  the Matsuno backward  time differencing           
c  scheme with semi-lagrangian advection is used here.           
c  the  semi-lagrangian  scheme   follows  the  method           
c  proposed by  Krishnamurti (1969) and mathur (1970).           
c  The  time  integration consists of a trial  forward           
c  step  and  a  corrector  step.  a  correction  term           
c  following   kanamitsu  (1974)  is  applied  to  the           
c  forcing functions during the corrector stage of the           
c  computation.                                                  
c                                                                               
      parameter (lg = 114,mg = 17)                                                    
c                                                                               
      common /a01/  l,      m,       l1,      m1,       m2,      g,             
     1              slat,   dphi                                                
      common /a02/  jin,    jout,    itape,   jb,       jc,      jd             
      common /a03/  ijk,    itcnt,   itime,   nstep,   iterr                    
      common /a04/  time,   dt,      nohist,  ioutpt,   nsbd,    alfa           
      common /a05/  tfcst,  toint,   niter,   nout,     title(3)                
      common /a06/  denom(6,mg), dx(mg),       cor(mg),    dy                   
      common /a07/  dx2iv(mg),   dxsqiv(mg),   dy2iv,      dysqiv               
      common /a08/  zm(lg,mg),   dhdx(lg,mg),  dhdy(lg,mg)                      
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg)                         
      common /a10/  up(lg,mg),   vp(lg,mg),    zp(lg,mg)                        
      common /a11/  uo(lg,mg),   vo(lg,mg),    zo(lg,mg)                        
      common /a12/  a(lg,mg),    b(lg,mg),     c(lg,mg)                         
      dimension  eu(lg,mg), ev(lg,mg), ez(lg,mg)                                
      equivalence (a,eu), (b,ev), (c,ez)
      integer debug
c                                                                               
      jinp      = jin+1                                                              
      m3        = m-3                                                                
      balfa     = (1.-alfa)*0.5                                                      
c                                                                               
c                                                                               
      do 11100  itime = jinp, jout                                                   
c                                                                             
c  Corrector  computation                 
c                                                                            
c                                                                               
c  Define some parameters required  for  iterative              
c  iterative initialization .           
c                                                                               
      if ( ijk.eq.0 ) then                                                       
      iii       = (itime-1)/nstep+1                                              
      iii       = mod(iii,2)+1                                                   
      dt        = abs(dt)*(-1.)**iii                                             
      endif                                                                   
      gdt       = g*dt                                                             
      time      = time + dt                                                        
c                                                                               
      if (nsbd.ne.0) then                                                      
      if (ijk.ne.0 ) then                                                     
      do 11102 j= 2, m1, m1-2                                                  
      do 11102 i= 1, l                                                      
         u(i,j) = alfa*u(i,j)+balfa*(u(i,j-1)+u(i,j+1))                 
         v(i,j) = alfa*v(i,j)+balfa*(v(i,j-1)+v(i,j+1))                 
11102    z(i,j) = alfa*z(i,j)+balfa*(z(i,j-1)+z(i,j+1))                 
      endif                                                                 
      endif                                                                   
c                                                                               
c  The forcing at the  southern  and  northern              
c  boundaries  at  time  t,  are  set  to  0.0                
c                                                                               
      do 11104 j= 1, m, m1                                                          
      do 11104 i= 1, l                                                           
         a(i,j) = 0.                                                        
         b(i,j) = 0.                                                        
11104    c(i,j) = 0.                                                        
c                                                                               
c  Computation of the forcing functions                
c  fort he other grid points.             
c                                                                               
      dyy       = dy2iv                                                               
      do 11106 j= 2, m1                                                            
         jp1    = j+1                                                             
         jm1    = j-1                                                             
         dxx    = dx2iv(j)                                                       
      do 11106 i= 1, l                                                           
         ip1    = i+1                                                         
         im1    = i-1                                                         
         if (ip1.gt.l) ip1 = 2                                                 
         if (im1.lt.1) im1 = l1                                                
         a(i,j) = -g*(z(ip1,j)-z(im1,j))*dxx + cor(j)*v(i,j)                
     &             -g*dhdx(i,j)                                              
         b(i,j) = -g*(z(i,jp1)-z(i,jm1))*dyy - cor(j)*u(i,j)                
     &             -g*dhdy(i,j)                                              
         c(i,j) = -z(i,j)*((u(ip1,j)-u(im1,j))*dxx                          
     &             +(v(i,jp1)-v(i,jm1))*dyy)                                 
11106 continue                                                                
c                                                                               
      call INTERP2                                                              
c                                                                               
c  Corrector step 
c                                                                               
c  First  the  correction terms to  be added  to the               
c  forcing functions are computed.                                 
c                                                                               
      dyy       = dysqiv                                                              
      do 11108 j= 3, m-2                                                           
         dxx    = dxsqiv(j)                                                      
         jp1    = j+1                                                            
         jp2    = j+2                                                             
         jm1    = j-1                                                            
         jm2    = j-2                                                            
      do 11108 i= 1, l                                                          
         zzz    = z(i,j)                                                      
         ip1    = i+1                                                         
         im1    = i-1                                                         
         if (ip1.gt.l) ip1 = 2                                                 
         if (im1.lt.1) im1 = l1                                                
         ip2    = i+2                                                         
         im2    = i-2                                                         
         if (ip2.gt.l) ip2 = i+3-l                                             
         if (im2.lt.1) im2 = l-3+i                                             
         eu(i,j)=      gdt*zzz*dxx*(-0.25*(u(ip2,j) + u(im2,j)                 
     &           -2.*u(i,j)) + u(ip1,j) + u(im1,j) -2.*u(i,j))           
         ev(i,j)=      gdt*zzz*dyy*(-0.25*(v(i,jp2) + v(i,jm2)                 
     &           -2.*v(i,j)) + v(i,jp1) + v(i,jm1)- 2.*v(i,j))           
         ez(i,j)= gdt*(dxx*(0.25*(z(ip2,j) + z(im2,j) -2.*zzz)            
     &            -(z(ip1,j) + z(im1,j) - 2.*zzz))+dyy*(0.25*(
     &            z(i,jp2) + z(i,jm2) - 2.*zzz)-(z(i,jp1) + z(
     &            i,jm1) - 2.*zzz)))                       
11108 continue                                                               
c                                                                               
      dxx       = dxsqiv(2)                                                       
      dyy       = dysqiv                                                          
      do 11110 i= 1, l                                                           
         zzz    = z(i,2)                                                       
         ip1    = i+1                                                          
         im1    = i-1                                                          
         if (ip1.gt.l) ip1 = 2                                                  
         if (im1.lt.1) im1 = l1                                                 
         ip2    = i+2                                                          
         im2    = i-2                                                          
         if (ip2.gt.l) ip2 = i+3-l                                              
         if (im2.lt.1) im2 = l-3+i                                              
         eu(i,2)= gdt*zzz*dxx*(-0.25*(u(ip2,2) + u(im2,2)-2.*
     &             u(i,2)) + u(im1,2) + u(ip1,2) - 2.*u(i,2))           
         ev(i,2)=  gdt*zzz*dyy*(v(i,3) + v(i,1) - 2.*v(i,2) +
     &        0.5*(v(i,2) - v(i,1)) - 0.25*(v(i,4) - v(i,2)))         
         ez(i,2)= gdt*(dxx*(0.25*(z(ip2,2)+z(im2,2) - 2.*zzz)            
     &            -(z(ip1,2) + z(im1,2) - 2.*zzz))+dyy*(0.25*
     &            (z(i,4) - zzz) - 0.5*(zzz - z(i,1))-(z(i,3)
     &            + z(i,1) - 2.*zzz)))                            
11110 continue                                                                
c                                                                               
      dxx       = dxsqiv(m1)                                                       
      dyy       = dysqiv                                                           
      do 11112 i= 1, l                                                            
           zzz  = z(i,m1)                                                       
           ip1  = i+1                                                           
           im1  = i-1                                                           
           if (ip1.gt.l) ip1 = 2                                                   
           if (im1.lt.1) im1 = l1                                                  
           ip2  = i+2                                                           
           im2  = i-2                                                           
           if (ip2.gt.l) ip2 = i+3-l                                               
           if (im2.lt.1) im2 = l-3+i                                               
      eu(i,m1)  = gdt*zzz*dxx*(-0.25*(u(ip2,m1)+u(im2,m1)                  
     &       -2.*u(i,m1))+u(ip1,m1)+u(im1,m1)-2.*u(i,m1))                
      ev(i,m1)  =  gdt*zzz*dyy*(v(i,m)+v(i,m2)-2.*v(i,m1)                   
     &      -0.5*(v(i,m)-v(i,m1))+0.25*(v(i,m1)-v(i,m3)))            
      ez(i,m1)  =  gdt*(dxx*(0.25*(z(ip2,m1)+z(im2,m1)-2.
     &       *zzz)-(z(ip1,m1)+z(im1,m1)-2.*zzz))+dyy*(-0.
     &       25*(zzz-z(i,m3))+0.5*(z(i,m)-zzz)-(z(i,m)+z(
     &       i,m2)-2.*zzz)))                               
11112 continue                                                                
c                                                                               
c  Start the corrector loop computation .             
c                                                                               
      dyy       = dy2iv                                                               
      do 11114 j= 2, m1                                                           
         dxx    = dx2iv(j)                                                      
         jp1    = j+1                                                           
         jm1    = j-1                                                           
      do 11114 i= 1, l                                                         
         ip1    = i+1                                                         
         im1    = i-1                                                         
         if (ip1.gt.l) ip1 = 2                                                 
         if (im1.lt.1) im1 = l1                                                
         fuq    = -g*(zo(ip1,j) - zo(im1,j))*dxx +
     &        cor(j)*vo(i,j)+eu(i,j) - g*dhdx(i,j)                                     
         fvq    = -g*(zo(i,jp1) - zo(i,jm1))*dyy -
     &        cor(j)*uo(i,j)+ev(i,j) - g*dhdy(i,j)                                     
         fzq    = -zo(i,j)*((uo(ip1,j)-uo(im1,j))*
     &      dxx+(vo(i,jp1)-vo(i,jm1))*dyy+ez(i,j))
c                        
         u(i,j) = up(i,j) + fuq*dt                                         
         v(i,j) = vp(i,j) + fvq*dt                                         
         z(i,j) = zp(i,j) + fzq*dt                                         
11114 continue                                                                
c                                                                               
c  After each time step, the day,  hour,  and  minute                
c  of the time step computation completed is  written              
c  into the listing file.This is optional if iprint=1            
c
      iprint    = 1
      if (debug.ge.1) then                                                                                     
      iday      = ifix (time/86400.)                                                
      ihour     = ifix ((time-iday*86400.)/3600.)                                   
      imin      = ifix ((time-iday*86400.-ihour*3600.)/60.)                         
      write (6,1000) iday, ihour, imin                                   
 1000 format(10x,' day =',i5,2x,' hour =',i5,2x,' min =',i5,)                  
      endif                                                 
11100 continue                                                                  
      return                                                                    
      end
c                                                                       
      function FINTP (f,i)                                                       
c                                                                               
c  This function computes the values of the  variables           
c  at  point  'p'   at  time  't'  using  the  9-point           
c  Lagrange interpolation scheme.                                
c                                                                               
      parameter (lg = 114,mg= 17)                                                    
c                                                                               
      common /a01/ l, m, l1, m1, m2, g, slat,dphi             
      common /a02/ jin, jout, itape, jb, jc, jd             
      common /gg1/ g1, g2, g3, g4, g5, g6, g7, g8, g9                          
      dimension    f(lg,mg)                                                    
c                                                                               
      ip1        = i+1                                                                 
      im1        = i-1                                                                 
      if (ip1.gt.l) ip1 = 2                                                        
      if (im1.lt.1) im1 = l1                                                       
      fintp      = g1*f(im1,jb)+g2*f(i,jb)+g3*f(ip1,jb)+
     &             g4*f(im1,jc)+g5*f(i,jc)+g6*f(ip1,jc)+
     &             g7*f(im1,jd)+g8*f(i,jd)+g9*f(ip1,jd)                                                      
      return                                                                    
      end
c                                                                       
      subroutine INDATA                                                          
c                                                                               
c  This routine  reads in the  input variables u, v, z,          
c  h, dh/dx, dh/dy.           
c                                                                               
      parameter ( lg = 114, mg = 17)                                                    
c                                                                               
      common /a01/  l,m,    l1,m1,  m2,g,  slat,dphi                                                
      common /a02/  jin,    jout,    itp,     jb,       jc,      jd             
      common /a03/  ijk,    itcnt,   itime,   nstep,   iterr                    
      common /a04/  time,   dt,      nohist,  ioutpt,   nsbd,    alfa           
      common /a08/  zm(lg,mg),   dhdx(lg,mg),  dhdy(lg,mg)                      
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg)                         
c                                                                               
      if (iterr.ne.0 ) then                                                       
        read (21,878) ((zm(i,j),i=1,l),j=1,m)                                    
        read (21,878) ((dhdx(i,j),i=1,l),j=1,m)                                  
        read (21,878) ((dhdy(i,j),i=1,l),j=1,m)                                  
  878   format(6e13.6)                                                          
      else                                                                      
      do 11120 i= 1, l                                                             
      do 11120 j= 1, m                                                           
      zm(i,j)   = 0.                                                       
      dhdx(i,j) = 0.                                                     
11120 dhdy(i,j) = 0.                                                     
      endif                                                                     
      itp       = 22                                                                    
      do 11121 kk= 1, itcnt                                                         
      read (itp,878) ((u(i,j),i=1,l),j=1,m)                                    
      read (itp,878) ((v(i,j),i=1,l),j=1,m)                                    
      read (itp,878) ((z(i,j),i=1,l),j=1,m)                                    
11121 continue                                                                  
      return                                                                    
      end
c                                                                       
      subroutine INTERP2                                                          
c                                                                               
c  This   subroutine   performs   the   semi-lagrangian          
c  advection computation.  the  method used is  similar          
c  to that proposed by Krishmamurti (1969)  and  mathur          
c  (1970).  the  9-point Lagrange interpolation  scheme          
c  is used for obtaining the  values at  the  departure          
c  points                                                        
c                                                                               
      parameter (lg = 114,mg = 17)                                                    
c                                                                               
      common /a01/  l,m,    l1,m1,   m2,g,  slat,dphi                                                
      common /a02/  jin,    jout,    itape,   jb,       jc,      jd             
      common /a04/  time,   dt,      nohist,  ioutpt,   nsbd,    alfa           
      common /a06/  denom(6,mg), dx(mg),       cor(mg),    dy                   
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg)                         
      common /a10/  up(lg,mg),   vp(lg,mg),    zp(lg,mg)                        
      common /a11/  uo(lg,mg),   vo(lg,mg),    zo(lg,mg)                        
      common /a12/  a(lg,mg),    b(lg,mg),     c(lg,mg)                         
      common /gg1/  g1, g2, g3, g4, g5, g6, g7, g8, g9                          
c                                                                               
      do 11130 j= 2, m1                                                              
         jb     = j-1                                                                
         jc     = j                                                                  
         jd     = j+1                                                                
         jp1    = j+1                                                                
         jm1    = j-1                                                                
      do 11131 i= 1, l                                                            
      do 11132 kk= 1, 2                                                         
         if (kk.eq.1) then                                                     
         uu     = u(i,jc)                                                       
         vv     = v(i,jc)                                                       
         au     = a(i,jc)                                                       
         bv     = b(i,jc)                                                       
         endif                                                                
         xp     = -(uu+.5*au*dt)*dt                                               
         yp     = -(vv+.5*bv*dt)*dt                                               
         yk1    = yp*(yp-dy)                                                     
         yk2    = (yp+dy)*(yp-dy)                                                
         yk3    = (yp+dy)*yp                                                     
         g1     = xp*(xp-dx(jm1))*yk1*denom(1,j)                                 
         g2     = (xp+dx(jm1))*(xp-dx(jm1))*yk1*denom(2,j)                       
         g3     = (xp+dx(jm1))*(xp*yk1*denom(1,j))                               
         g4     = xp*(xp-dx(j))*yk2*denom(3,j)                                   
         g5     = (xp+dx(j))*(xp-dx(j))*yk2*denom(4,j)                           
         g6     = (xp+dx(j))*xp*yk2*denom(3,j)                                   
         g7     = (xp-dx(jp1))*xp*yk3*denom(5,j)                                 
         g8     = (xp+dx(jp1))*(xp-dx(jp1))*yk3*denom(6,j)                       
         g9     = (xp+dx(jp1))*(xp*yk3*denom(5,j))                               
c                                                                               
c  Computations of the values of u, v, a, b at point p                
c  are  accomplished  by calling the function 'FINTP.'              
c                                                                               
         ipass  = i                                                           
         xx     = fintp(u(1,1),ipass)                                            
         xy     = fintp(v(1,1),ipass)                                            
         ax     = fintp(a(1,1),ipass)                                            
         bx     = fintp(b(1,1),ipass)                                            
c                                                                               
c  Redefine u,v,a,b at p                                      
c                                                                               
         if (kk.eq.1) then                                                    
         uu     = xx                                                           
         vv     = xy                                                           
         au     = ax                                                           
         bv     = bx                                                           
         else                                                                
         up(i,j)= xx                                                      
         vp(i,j)= xy                                                      
c                                                                               
c  Calculation of z value at point p.                               
c                                                                               
          xz    = fintp(z(1,1),ipass)                                          
          zp(i,j)= xz                                                      
          cx    = fintp(c(1,1),ipass)                                          
          uo(i,j)= xx + ax*dt                                              
          vo(i,j)= xy + bx*dt                                              
          zo(i,j)= xz + cx*dt                                              
          endif                                                               
11132 continue                                                             
11131 continue                                                               
11130 continue                                                                  
      return                                                                    
      end
c                                                                       
      subroutine INVART(debug)                                                          
c                                                                               
c  This subroutine computes  the closed  domain invariants          
c  of the single level primitive equation model.  these          
c  are the domain average (1) geopotential height,  (2)          
c  potential vorticity, (3) potential vorticity square,          
c  and (4) ln(z**2/2) .the root mean square  divergence          
c  is also computed.                                             
c                                                                               
      parameter (lg = 114,mg = 17)                                                    
c                                                                               
      common /a01/  l,m,    l1,m1,   m2,g,  slat,dphi                                                
      common /a02/  jin,    jout,    itape,   jb,       jc,      jd             
      common /a03/  ijk,    itcnt,   itime,   nstep,    iterr                   
      common /a06/  denom(6,mg), dx(mg),       cor(mg),    dy                   
      common /a07/  dx2iv(mg),   dxsqiv(mg),   dy2iv,      dysqiv               
      common /a08/  zm(lg,mg),   dhdx(lg,mg),  dhdy(lg,mg)                      
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg) 
      integer debug
c                                                                               
      areaiv    = 1./float(l1*m2)                                                  
      sum       = 0.                                                                  
      sumz      = 0.                                                                  
      sumzsq    = 0.                                                               
      sume      = 0.                                                                
      pvsq      = 0.                                                                
      divsq     = 0.                                                                
c                                                                               
      do 11140 j= 2, m1                                                              
      do 11140 i= 1, l1                                                            
         ip1    = i+1                                                             
         im1    = i-1                                                             
         if (im1.lt.1) im1 = l1                                                   
         potvor = ((v(ip1,j)-v(im1,j))/(2.*dx(j))-(u(i,j+1)-
     &                      u(i,j-1))/(2.*dy)+cor(j))/z(i,j)                  
         pvsq   = pvsq + potvor*potvor                                         
         sum    = sum  + potvor                                                 
         sume   = sume + z(i,j)*0.5*(u(i,j)*u(i,j)+v(i,j)*v(i,j)               
     &          + g*z(i,j))                                                  
         if (iterr.eq.1) sume = sume+z(i,j)*g*zm(i,j)                              
         sumzsq = sumzsq + alog(z(i,j)*z(i,j)/2)                                
         sumz   = sumz   + z(i,j)                                                  
         div    = (u(ip1,j)-u(im1,j))*dx2iv(j)+(v(i,j+1)-v(i,j-1)
     &             )*dy2iv                                     
         divsq  = divsq + div*div                                               
11140 continue                                                                  
c                                                                               
      sum       = sum * areaiv                                                     
      divsq     = sqrt(divsq*areaiv)                                               
      sume      = sume * areaiv                                                    
      sumz      = sumz * areaiv                                                    
      sumzsq    = sumzsq*areaiv                                                    
      pvsq      = pvsq * areaiv                                                    
      if (debug.ge.1) then
      write (6,50) sumz, sum, pvsq, sumzsq, sume, divsq                          
  50  format(//,1x,'mean geopotential  height  ',e13.4,/,                       
     &          1x,'mean pot. vorticity        ',e13.4,/,                       
     &          1x,'mean square pot. vorticity ',e13.4,/,                       
     &          1x,'mean ln(0.5*z*z)           ',e13.4,/,                       
     &          1x,'mean total energy          ',e13.4,/,                       
     &          1x,'root mean square divergence',e13.4,/)
      endif
      return                                                                    
      end                                                                       
c
      subroutine ZEROS                                                           
c                                                                               
c  This subroutine sets the values of the variables
c  in the common blocks to zeros .                                  
c                                                                               
      parameter (lg = 114,mg = 17)                                                    
c                                                                               
      common /a05/  tfcst,  toint,   niter,   nout,     title(3)                
      common /a06/  denom(6,mg), dx(mg),       cor(mg),    dy                   
      common /a07/  dx2iv(mg),   dxsqiv(mg),   dy2iv,      dysqiv               
      common /a08/  zm(lg,mg),   dhdx(lg,mg),  dhdy(lg,mg)                      
      common /a09/  u(lg,mg),    v(lg,mg),     z(lg,mg)                         
      common /a10/  up(lg,mg),   vp(lg,mg),    zp(lg,mg)                        
      common /a11/  uo(lg,mg),   vo(lg,mg),    zo(lg,mg)                        
      common /a12/  a(lg,mg),    b(lg,mg),     c(lg,mg)                         
c                                                                               
      do 11150 i= 1, 3                                                                
      title(i)  = 0.                                                           
11150 continue                                                                  
      do 11151 i= 1, 6                                                               
      do 11151 j= 1, mg                                                            
11151 denom(i,j)= 0.                                                       
      do 11152 i= 1, mg                                                              
      dx(i)     = 0.                                                           
      cor(i)    = 0.                                                           
      dx2iv(i)  = 0.                                                          
      dxsqiv(i) = 0.                                                          
      do 11152 j= 1, lg                                                         
      zm(j,i)   = 0.                                                       
      dhdx(j,i) = 0.                                                       
      dhdy(j,i) = 0.                                                       
      u(j,i)    = 0.                                                         
      v(j,i)    = 0.                                                         
      z(j,i)    = 0.                                                         
      up(j,i)   = 0.                                                         
      vp(j,i)   = 0.                                                         
      zp(j,i)   = 0.                                                         
      uo(j,i)   = 0.                                                         
      vo(j,i)   = 0.                                                         
      zo(j,i)   = 0.                                                         
      a(j,i)    = 0.                                                          
      b(j,i)    = 0.                                                          
      c(j,i)    = 0.                                                          
11152 continue                                                                  
      return                                                                    
      end                                                                       
