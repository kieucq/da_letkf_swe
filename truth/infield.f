        program infield                                                         
c                                                                               
c  this  program  can  be  used  either  to prepare the                         
c  initial streamfunction  field for  the non-divergent                         
c  barotropic  model  forecast  or the u,v,z fields for                         
c  the shallow water model iterativeinitialization/forecast.                    
c  if the program is to  be used for preparing  the initial                     
c  field  for  the  barotropic model, set parameters                            
c    iopt1  =  0, and                                                           
c    iopt2  =  0                                                                
c  in this case the program reads in the  initial wind                          
c  field  from unit 21 and outputs  the streamfunction                          
c  field to unit 22.                                                            
c  if the program is to be used for preparing the initial                       
c  field for the shallow water model without terrain,set                        
c    iopt1  =  1, and                                                           
c    iopt2  =  0                                                                
c  the  program will read in the initial u,v field  from                        
c  unit 21 and outputs the u,v,z field to unit 22.                              
c  for the preparation of initial field  for shallow water                      
c  model with terrain, set                                                      
c    iopt1   = 1 ,and                                                           
c    iopt2   = 1 .                                                              
c                                                                               
c  originator :       the florida state university. last revison by             
c                     L. bounoua , july 1994.                                   
c                                                                               
c  input      : (1) initial   grid   point   zonal   and  meridional            
c                   components of the wind fields from unit 21.                 
c               (2) if  program is  used for  preparing the  initial            
c                   fields  for the shallow water model with terrain,           
c                   then the  reduced  terrain  field and  its zonal            
c                   and meridional gradients  are read in from  unit            
c                   20.                                                         
c                                                                               
c  output     : (1) logical  unit  22  contains the   output  fields            
c                   described above.                                            
c               (2) logical unit 10 is the dayfile file.                        
c               (3) logical unit 23 contains terrain,zonal and its mer-         
c                   idional gradients .                                         
c                                                                               
c  subroutines :    bound, const,cycle, jacmod ,lapmod , relax, zfield, terr.   
c  called           terr ,stream,zfield .                                       
c                                                                               
c  definition  :                                                                
c                                                                               
c  l           : east-west   dimension of the domain                            
c  m           : north-south dimension of the domain                            
c  ladd        : number of points added to make the domain cyclic               
c  u,v         : horizontal wind components                                     
c  psi         : streamfunction                                                 
c  z           : height field                                                   
c  cor         : coriolis parameter for each latitude                           
c  dx          : east-west grid spacing for each latitude                       
c  dx2         : east-west grid spacing for terrain field                       
c  aaa,uu,vv   : work arrays                                                    
c  wrka1,wrka  : work arrays                                                    
c  hd,dhdy,dhdx: terrain and its gradients                                      
c                                                                               
c                                                                               
      parameter(l = 108,m = 17,ladd = 6,lcy=l+ladd,iopt1 = 1,iopt2 = 0 )        
      real    u(l,m)     ,  v(l,m)      , psi(l,m)                              
      real    aaa(lcy,m) ,  uu(lcy,m)   , vv(lcy,m)                             
      real    z(lcy,m)   ,  wrka1(lcy,m), wrka(lcy,m)                           
      real    cor(m)     ,  dx(m)       , dx2(m)                                
      real    dhdy(lcy,m),  dhdx(lcy,m) , hd(lcy,m)                             
      real    az(l,m)                                                           
      equivalence (wrka,vv)                                                     
c                                                                               
c  open the input-output files .unit 20 and 23 are opened only                  
c  if the terrain is used .                                                     
c                                                                               
      open (21,file='infield_in.dat',status='old')                                   
      open (22,file='infield_out.dat',status='unknown')                         
      if   (iopt2.ne.0) then                                                    
      open (20,file='terrain.dat',status='old')                                 
      open (23,file='terrain.out',status='unknown')                             
      endif                                                                     
c                                                                               
c  print model functions                                                        
c                                                                               
      print *,' '                                                               
      print *,'    model status'                                                
      print *,' '                                                               
      if (iopt1.eq.0.and.iopt2.eq.0) then                                       
      print *,'    initialization for barotropic model '                        
      print *,' '                                                               
      print *,'    the model reads the initial u,v field from unit 21'          
      print *,'    and outputs the streamfunction field to unit 22'             
      print *,' '                                                               
      else                                                                      
      if (iopt1.eq.1.and.iopt2.eq.0) then                                       
      print *,'    initialization for the single level '                        
      print *,'    primtive equation model without terrain'                     
      print *,' '                                                               
      print *,'    the model reads in the input u,v field from unit 21'         
      print *,'    and outputs the u,v,z fields to unit 22'                     
      print *,' '                                                               
      else                                                                      
      if (iopt1.eq.1.and.iopt2.eq.1) then                                       
      print *,'    initialization for the single level '                        
      print *,'    primtive equation model with terrain'                        
      print *,' '                                                               
      print *,'    the model reads in the input u,v field from unit 21'         
      print *,'    and outputs the u,v,z fields to unit 22'                     
      print *,'    the initial terrain field is read from unit 20'              
      print *,'    and the output field is written to unit 23'                  
      print *,' '                                                               
      endif                                                                     
      endif                                                                     
      endif                                                                     
      print *,' '                                                               
      print *,'                        note  : '                                
      print *,'    the grid size and/or the southernmost latitude need'         
      print *,'    to be changed in subroutine "const" if the         '         
      print *,'    domain/or data is changed'                                   
      print *,' '                                                               
      print *,' '                                                               
c                                                                               
c  rewind all input tapes .                                                     
c                                                                               
      if(iopt2.ne.0)rewind (20)                                                 
      rewind (21)                                                               
      rewind (22)                                                               
c                                                                               
      l1        = l-1                                                           
      m1        = m-1                                                           
      m2        = m-2                                                           
c                                                                               
c  call subroutine  const  to  define  the constants                            
c  required by the program .                                                    
c                                                                               
      call const (m,dx,dy,cor,beta,slat,dphi)                                   
c                                                                               
c     read in the initial u, v fields .                                         
c                                                                               
      print *,'printing first and last value of each field  '                   
      read (21,878) ((az(ii,jj),ii=1,l),jj=1,m)                                 
      read (21,878) ((u (ii,jj),ii=1,l),jj=1,m)                                 
      read (21,878) ((v (ii,jj),ii=1,l),jj=1,m)                                 
  878 format(6e13.6)                                                            
c                                                                               
c  print first and last elements of u and v arrays as a check.                  
c                                                                               
      write (6,1000) u(1,1),u(l,m)                                              
      write (6,1001) v(1,1),v(l,m)                                              
      write (6,1002)                                                            
 1000 format(1x,'u(1,1)=',e12.4,5x,'u(l,m)=',e12.4)                             
 1001 format(1x,'v(1,1)=',e12.4,5x,'v(l,m)=',e12.4)                             
 1002 format(//)                                                                
 1003 format(1x,'z(1,1)=',e12.4,5x,'z(l,m)=',e12.4)                             
c                                                                               
c  call subroutine  streamf to compute  the streamf-                            
c  unction field.                                                               
c                                                                               
                                                                                
      call STREAMF (u,v,l,lcy,l,m,dy,dx,wrka,psi)                               
                                                                                
c                                                                               
c   call subroutine  cycle to  puts east-west cyclic                            
c   boundary condition on the streamfunction field .                            
c                                                                               
      do 10100 j = 1, m                                                         
      do 10100 i = 1, l                                                         
         aaa(i,j)= psi(i,j)                                                     
10100 continue                                                                  
      call cycle (aaa,lcy,m,1)                                                  
      if (iopt1.eq.0) then                                                      
      write (22,879) ((aaa(i,j),i=1,lcy),j=1,m)                                 
  879 format(6e13.6)                                                            
      else                                                                      
c                                                                               
c  call subroutine  zfield to  compute the height field                         
c  from the streamfunction field.                                               
c  the mean height filed of the 700 mb surface is considered .                  
c                                                                               
      zbar      = 3160.                                                         
      zbar2     = 0.                                                            
      call zfieldmod (aaa,uu,vv,lcy,m,dx,dy,cor,wrka1,z,1.,1.,zbar,beta)        
      do 10102 j= 1, m                                                          
      do 10102 i= 1, lcy-1                                                      
         zbar2  = z(i,j) + zbar2                                                
10102 continue                                                                  
      zbar2     = zbar2/(m*(lcy-1))                                             
      do 10104 j= 1, m                                                          
      do 10104 i= 1, lcy                                                        
         z(i,j) = z(i,j)-zbar2+zbar                                             
10104 continue                                                                  
c                                                                               
c  call subroutine cycle to put east - west cyclic conditions on                
c  the u,v fields .                                                             
c                                                                               
      do 10105 j= 1, m                                                          
      do 10105 i= 1, lcy                                                        
        aaa(i,j)= 0.                                                            
10105 continue                                                                  
      do 10106 j= 1, m                                                          
      do 10106 i= 1, l                                                          
        aaa(i,j)= u(i,j)                                                        
10106 continue                                                                  
                                                                                
      call cycle (aaa,lcy,m,1)                                                  
      write (22,879) ((aaa(i,j),i=1,lcy),j=1,m)                                 
                                                                                
      do 10107 j= 1, m                                                          
      do 10107 i= 1, lcy                                                        
        aaa(i,j)= 0.                                                            
10107 continue                                                                  
      do 10108 j= 1, m                                                          
      do 10108 i= 1, l                                                          
      aaa(i,j)  = v(i,j)                                                        
10108 continue                                                                  
      call cycle (aaa,lcy,m,1)                                                  
                                                                                
      write (22,879) ((aaa(i,j),i=1,lcy),j=1,m)                                 
c                                                                               
c  read the terrain field if it is used. read in terrain                        
c  field at l by m and then run cycle to get cyclic values                      
c  between l and lcy in zonal direction                                         
c                                                                               
      if (iopt2.eq.1) then                                                      
      read (20,880) ((wrka1(ii,jj),ii=1,l),jj=1,m)                              
      call cycle (wrka1,lcy,m,1)                                                
      aamax     = wrka1(1,1)                                                    
c                                                                               
c  compute the maximum valueof the terrain field and                            
c  normalized it .                                                              
c                                                                               
      do 10110 j= 1, m                                                          
      do 10110 i= 1, lcy                                                        
10110 if (wrka1(i,j).gt.aamax) aamax = wrka1(i,j)                               
      do 10112 j= 1, m                                                          
      do 10112 i= 1, lcy                                                        
10112 wrka1(i,j)= (wrka1(i,j)/aamax) * 1000.0                                   
c                                                                               
c  subroutine terr create tape21 for shallow water                              
c  model when iopt2.ne.0                                                        
c                                                                               
      call terr (wrka1,lcy,m,dphi,slat,dhdy,dhdx,dx2,hd)                        
  880 format(6e13.6)                                                            
      do 10114 j= 1, m                                                          
      do 10114 i= 1, lcy                                                        
10114    z(i,j) = z(i,j)-1.0*wrka1(i,j)                                         
      endif                                                                     
c                                                                               
c  write output fields                                                          
c                                                                               
      write (22,879) ((z(i,j),i=1,lcy),j=1,m)                                   
      write (6,1000) u(1,1),u(l,m)                                              
      write (6,1001) v(1,1),v(l,m)                                              
      write (6,1003) z(1,1),z(l,m)                                              
      write (6,1002)                                                            
      endif                                                                     
      stop                                                                      
      end                                                                       
c                                                                               
      subroutine bound (a,l,m,n,l1,m1,m2,ladd)                                  
c                                                                               
c  this  subroutine  is  used to  obtain  the southern                          
c  and northern boundary values by extrapolation.  for                          
c  domain  that is  cyclic in the zonal direction, the                          
c  values  are  extrapolated  linearly  outwards . for                          
c  domain  that  have  zonal boundary values formed by                          
c  extending  the  values  1  grid point  outwards,  a                          
c  similar extension is done here for the southern and                          
c  northern boundary values.                                                    
c                                                                               
c  definitions  :                                                               
c    a(l,m,n)   : variable  whose  southern  and  northern                      
c                 boundary  values have to be  obtained by                      
c                 extrapolation.                                                
c    l          : number  of  grid   points  in  the zonal                      
c                 direction.                                                    
c    m          : number  of grid points in the meridional                      
c                 direction.                                                    
c    n          : number of levels in the vertical.                             
c    l1         : l-1                                                           
c    m1         : m-1                                                           
c    m2         : m-2                                                           
c    ladd       : number  of  grid  points  added  to  the                      
c                 eastern  boundary to   make  the  domain                      
c                 cyclic  in  the zonal direction.                              
c                                                                               
      real a(l,m,n)                                                             
      k         = 1                                                             
c                                                                               
      do 10110 i= 1, l                                                          
      a(i,1,k)  = 2.*a(i,2,k)-a(i,3,k)                                          
      a(i,m,k)  = 2.*a(i,m1,k)-a(i,m2,k)                                        
10110 continue                                                                  
c                                                                               
c  the following  loop is for the non-cyclic case only.                         
c                                                                               
      if (ladd.eq.0) then                                                       
      do 10111 i= 1, l                                                          
      a(i,1,k)  = a(i,2,k)                                                      
10111 a(i,m,k)  = a(l1,m1,k)                                                    
      endif                                                                     
      return                                                                    
      end                                                                       
c                                                                               
      subroutine const (m,dx,dy,cor,beta,slat,dphi)                             
c                                                                               
c  this   subroutine   defines   the   constants   and                          
c  parameters required by program  infield . users  of                          
c  this program need only to modify the parameter                               
c  statement in program  infield  and this subroutine                           
c  to their requirements.                                                       
c                                                                               
c  definitions   :                                                              
c                                                                               
c  m             : number of grid points in the  meridional                     
c                  direction.                                                   
c  dx(m)         : grid  size in  meters  for  each  row of                     
c                  grid points in the zonal direction.                          
c  dy            : grid size in meters  for  the meridional                     
c                  direction.                                                   
c  cor(m)        : coriolis parameter for each  row of grid                     
c                  points.                                                      
c  beta          : beta plane parameter (=df/dy).                               
c                                                                               
c  omega         : earth rotation (rad/s)                                       
c  er            : earth radius   (m    )                                       
c                                                                               
      real cor(m),dx(m)                                                         
      dphi      = 2.5                                                           
      slat      = 0.0                                                           
      anlat     = slat+(m-1)*dphi                                               
      omega     = 7.292e-05                                                     
      er        = 6.37122e06                                                    
      pi        = 4.0*atan(1.0)                                                 
      rad       = pi/180.                                                       
      dy        = dphi*111.1*1000.0                                             
      phi       = slat                                                          
      do 10120 j= 1, m                                                          
         dx(j)  = dy*cos(phi*rad)                                               
         cor(j) = 2.0*omega*sin(phi*rad)                                        
         phi    = phi+dphi                                                      
10120 continue                                                                  
      the       = (slat+anlat)/2.0                                              
      beta      = 2.0*omega*cos(the*rad)/er                                     
      return                                                                    
      end                                                                       
c                                                                               
      subroutine cycle (z,l,m,n)                                                
c                                                                               
c  this  subroutine  is  used  for  interpolating  the                          
c  values over 6 grid points  at the eastern  boundary                          
c  to create a cyclic boundary in the zonal direction.                          
c                                                                               
c  definitions    :                                                             
c                                                                               
c   z(l,m,n)      : array  whose  zonal  values  are  to  be                    
c                   interpolated   over  6  grid  points  to                    
c                   create a cyclic boundary.                                   
c   l             : number  of  grid  points  in  the  zonal                    
c                   direction.                                                  
c   m             : number of grid points in the  meridional                    
c                   direction.                                                  
c   n             : number of levels in the vertical.                           
c                                                                               
      real z(l,m,n)                                                             
      k         = 1                                                             
      m1        = m-1                                                           
      l1        = l-1                                                           
      l2        = l-2                                                           
      l3        = l-3                                                           
      l4        = l-4                                                           
      l5        = l-5                                                           
      l6        = l-6                                                           
      l7        = l-7                                                           
      l8        = l-8                                                           
      do 10130 j= 1, m                                                          
      z(l,j,k)  = z(1,j,k)                                                      
      z(l3,j,k) = (z(l6,j,k)+z(1,j,k))/2.                                       
      z(l5,j,k) = (z(l3,j,k)+z(l6,j,k)*8.0-z(l7,j,k)*3.0)/6.                    
      z(l1,j,k) = (z(l3,j,k)+z(l,j,k)*8.0-z(2,j,k)*3.0)/6.                      
      z(l4,j,k) = (z(l6,j,k)*3.0+z(l3,j,k)*6.0-z(l1,j,k))/8.0                   
      z(l2,j,k) = (z(l1,j,k)*3.0+z(l3,j,k)*6.0-z(l5,j,k))/8.0                   
10130 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
      subroutine jacmod (a,b,c,dx,dy,l,m,n,l1,m1,l2,m2,ladd)                    
c                                                                               
c  this   subroutine   computes  the  advective term,                           
c  the    jacobian    of   (psi, delsquared psi + f),                           
c  following arakawa (1966).  a nine-point stencil is                           
c  used . note that this subroutine is similar in pri-                          
c  nciples to subroutine jac ,but the cyclicity condi-                          
c  tion has been added .                                                        
c                                                                               
c  defiitions      :                                                            
c                                                                               
c  a(l,m,n)        : this array  contains  the jacobian  of                     
c                    psi and delsquared psi + f, where f is                     
c                    the earth's  vorticity.                                    
c                                                                               
c  b(l,m,n)        : the streamfunction field,psi.                              
c                                                                               
c  c(l,m,n)        : the  absolute  vorticity or delsquared                     
c                    (psi + f).                                                 
c  dx(m)           : grid size in the  zonal direction  for                     
c                    each row of grid points in meters.                         
c  dy              : grid size in the  meridional direction                     
c                    in meters.                                                 
c  l               : number  of grid  points  in  the zonal                     
c                    direction.                                                 
c  m               : number   of   grid   points   in   the                     
c                    meridional direction.                                      
c  n               : number of levels in the vertical.                          
c  l1              : l-1                                                        
c  m1              : m-1                                                        
c  l2              : l-2                                                        
c  m2              : m-2                                                        
c  ladd            : number  of  grid  points  added to the                     
c                    eastern boundary  to  make  the domain                     
c                    cyclic in the zonal direction.                             
c                                                                               
      real a(l,m,n), b(l,m,n), c(l,m,n), dx(m)                                  
      k          = 1                                                            
      do 10140 j = 2, m1                                                        
         dm      = 12.*dx(j)*dy                                                 
         jp1     = j+1                                                          
         jm1     = j-1                                                          
      do 10140 i = 1, l                                                         
         if (i-1)80,80,81                                                       
   80    im1     = l1                                                           
         ip1     = 2                                                            
         go to 83                                                               
   81    if(i-l)82,80,80                                                        
   82    im1     = i-1                                                          
         ip1     = i+1                                                          
   83    continue                                                               
      a(i,j,k)   =(b(i,j-1,k)+b(ip1,j-1,k)-b(i,j+1,k)-b(ip1,j+1,k))             
     &  *(c(ip1,j,k)-c(i,j,k))+(b(im1,j-1,k)+b(i,j-1,k)-b(im1,j+1,k             
     &  )-b(i,j+1,k))*(c(i,j,k)-c(im1,j,k))+(b(ip1,j,k)+b(ip1,j+1,k             
     &  )-b(im1,j,k)-b(im1,j+1,k))*(c(i,j+1,k)-c(i,j,k))+(b(ip1,j-1             
     &  ,k)+b(ip1,j,k)-b(im1,j-1,k)-b(im1,j,k))*(c(i,j,k)-c(i,j-1,k             
     &  ))+(b(ip1,j,k)-b(i,j+1,k))*(c(ip1,j+1,k)-c(i,j,k))+(b(i,j-1             
     &  ,k)-b(im1,j,k))*(c(i,j,k)-c(im1,j-1,k))+(b(i,j+1,k)-b(im1,j             
     &  ,k))*(c(im1,j+1,k)-c(i,j,k))+(b(ip1,j,k)-b(i,j-1,k))*(c(i,j             
     &  ,k)-c(ip1,j-1,k))                                                       
c                                                                               
      a(i,j,k)  = a(i,j,k)/dm                                                   
10140 continue                                                                  
      do 10141 i= 1, l                                                          
      if (i-1)70,70,71                                                          
   70 im1       = l1                                                            
      ip1       = 2                                                             
      go to 73                                                                  
   71 if (i-l) 72,70,70                                                         
   72 im1       = i-1                                                           
      ip1       = i+1                                                           
   73   continue                                                                
c                                                                               
      a(i,1,k)  = (b(i,1,k)+b(ip1,1,k)-b(i,2,k)-b(ip1,2,k))*(c(i,1,k)+          
     &    c(ip1,1,k))-(b(im1,1,k)+b(i,1,k)-b(im1,2,k)-b(i,2,k))*(c(im1          
     &    ,1,k)+c(i,1,k))+(b(ip1,1,k)+b(ip1,2,k)-b(im1,1,k)-b(im1,2,k)          
     &    )*(c(i,1,k)+c(i,2,k))+(b(ip1,1,k)-b(i,2,k))*(c(i,1,k)+c(ip1,          
     &    2,k))+(b(i,2,k)-b(im1,1,k))*(c(im1,2,k)+c(i,1,k))                     
c                                                                               
      dm        = 12.*dx(1)*dy                                                  
      a(i,1,k)  = a(i,1,k)/dm                                                   
      a(i,m,k)  = (b(i,m-1,k)+b(ip1,m-1,k)-b(i,m,k)-b(ip1,m,k))*(c(i,m,         
     &    k)+c(ip1,m,k))-(b(im1,m-1,k)+b(i,m-1,k)-b(im1,m,k)-b(i,m,k))*         
     &    (c(im1,m,k)+c(i,m,k))-(b(ip1,m-1,k)+b(ip1,m,k)-b(im1,m-1,k)-b         
     &    (im1,m,k))*(c(i,m-1,k)+c(i,m,k))-(b(i,m-1,k)-b(im1,m,k))*(c(i         
     &    m1,m-1,k)+c(i,m,k))-(b(ip1,m,k)-b(i,m-1,k))*(c(i,m,k)+c(ip1,m         
     &    -1,k))                                                                
c                                                                               
      dm        = 12.*dx(m)*dy                                                  
      a(i,m,k)  = a(i,m,k)/dm                                                   
10141 continue                                                                  
c                                                                               
c  the next do loop is for the non-cyclic case only.                            
c                                                                               
      if (ladd.eq.0) then                                                       
      do 10142 j = 1, m                                                         
      a(1,j,k)   = a(2,j,k)                                                     
      a(l,j,k)   = a(l1,j,k)                                                    
10142 continue                                                                  
      endif                                                                     
      return                                                                    
      end                                                                       
c                                                                               
      subroutine lapmod(a,b,dx,dy,l,m,n,l1,m1,l2,m2,ladd)                       
c                                                                               
c  this subroutine lap calculates  the laplacian of any                         
c  scalar  field.  the  standard   second order  finite                         
c  differencing is used  .   two   forms   of  boundary                         
c  conditions are included,  one  that  is   cyclic  in                         
c  the zonal  direction,  and   the  other that is  non                         
c  cyclic  in  the  zonal  direction.  for  the  cyclic                         
c  case, the southern and northern  boundary values are                         
c  formed by linear extrapolation  outwards , while for                         
c  the non-cyclic case , the boundary values are formed                         
c  by extending the values one grid point outwards.                             
c    note that this subroutine is similar in principles                         
c  to subroutine jac ,but the cyclicity condition has be-                       
c  en added .                                                                   
c                                                                               
c  definitions    :                                                             
c                                                                               
c   a(l,m,n)      : laplacian  of  b.  for our  purpose this                    
c                   corresponds  to the  relative vorticity.                    
c   b(l,m,n)      : any scalar  field.  for our purpose this                    
c                   corresponds to the streamfunction field.                    
c                                                                               
c   dy            : grid   distance   in   meters   in   the                    
c                   meridional  direction.                                      
c   l             : number  of  grid  points  in  the  zonal                    
c                   direction.                                                  
c   m             : number of grid points in the  meridional                    
c                   direction.                                                  
c   n             : number of vertical levels.                                  
c   l1            : l-1                                                         
c   m1            : m-1                                                         
c   l2            : l-2                                                         
c   m2            : m-2                                                         
c   ladd          : number  of  grid  points   added to  the                    
c                   eastern  boundary  to  make  the  domain                    
c                   cyclic in the  zonal  direction.                            
c                                                                               
      real  a(l,m,n), b(l,m,n), dx(m)                                           
      k         = 1                                                             
      do 10150 j= 2, m1                                                         
         jp1    = j+1                                                           
         jm1    = j-1                                                           
      do 10150 i= 2, l1                                                         
10150 a(i,j,k)  = ((b(i+1,j,k)+b(i-1,j,k)-2.*b(i,j,k))/dx(j)**2                 
     &              +(b(i,jp1,k)+b(i,jm1,k)-2.*b(i,j,k))/dy**2)                 
      do 10151 j= 2, m1                                                         
      a(1,j,k)  =  ((b(2,j,k)+b(l1,j,k)-2.*b(1,j,k))/dx(j)**2                   
     &            +(b(1,j+1,k)+b(1,j-1,k)-2.*b(1,j,k))/dy**2)                   
10151 a(l,j,k)  = a(1,j,k)                                                      
c                                                                               
c  subroutine  bound  is  used   to set  the southern                           
c  and northern boundary conditions.                                            
c                                                                               
      call bound (a,l,m,n,l1,m1,m2,ladd)                                        
c                                                                               
c  the next do loop is for non-cyclic case only.                                
c                                                                               
      if (ladd.eq.0) then                                                       
      do 10152 j = 1, m                                                         
      a(1,j,k)   = a(2,j,k)                                                     
10152 a(l ,j,k)  = a(l1,j,k)                                                    
      endif                                                                     
      return                                                                    
      end                                                                       
c                                                                               
      subroutine relaxmod (x,zzinv,z,y,zinv,l,l1,m1,m,lcycle)                   
c                                                                               
c  this  subroutine  solves  the  poisson's  equation,                          
c  delsquared  x(i,j) = y(i,j),   using  a  successive                          
c  approximation method. the sequential overrelaxation                          
c  method is used. in program "infield", this  routine                          
c  is used to obtain the height filed from the reverse                          
c  balance laws.  except for some of the values of the                          
c  constants, this subroutine is similar to subroutine                          
c  relaxt  in program (baro).                                                   
c    note that this subroutine is similar in principles                         
c  to subroutine relax where the cyclicity condition has                        
c  been added .                                                                 
c                                                                               
c  definitions    :                                                             
c                                                                               
c  x(l,m)         : height field.                                               
c  zzinv          : 1/(dy**2),  where dy is the grid size in                    
c                   meters in the meridional direction.                         
c  z(m)           : dx(j)**2, where dx(j) is  the grid  size                    
c                   in meters for each row of grid points in                    
c                   the zonal direction.                                        
c  y(l,m)         : forcing function field.                                     
c  zinv(m)        : 1/z(j) where j=1...(1)...m.                                 
c  l              : number  of  grid  points  in  the  zonal                    
c                   direction.                                                  
c  l1             : l-1                                                         
c  m              : number  of grid points in the meridional                    
c                   direction.                                                  
c  m1             : m-1                                                         
c  lcycle         : control integer indicating  which domain                    
c                   in  the   zonal  direction  to   perform                    
c                   relaxation on.  if lcycle=1,  relaxation                    
c                   is  performed  for the domain, i=2...l1,                    
c                   j=2...m1;  if   lcycle=other   integers,                    
c                   relaxation  is performed for the domain,                    
c                   i=1...l,  j=2...m1.                                         
c                                                                               
      real  x(l,m),y(l,m),z(m),zinv(m)                                          
      npts      = l*(m-2)                                                       
      if (lcycle.eq.1) npts = (l-2)*(m-2)                                       
      mm        = 2                                                             
      mmm       = m1                                                            
      alfa      =.46                                                            
      ia        = 1000                                                          
      eps       = 1.                                                            
      nsc       = 0                                                             
      lsc       = -1                                                            
      lcyc1     = 1                                                             
      lcyc2     = l                                                             
      if (lcycle.eq.1) then                                                     
      lcyc1     = 2                                                             
      lcyc2     = l1                                                            
      endif                                                                     
      do 10160 j= 2, m1                                                         
      do 10160 i= 2, l1                                                         
10160    x(i,j) = 0.0                                                           
   15 nrel      = 0                                                             
      l1        = l-1                                                           
      do 10161 j= mm, mmm                                                       
         jp1    = j+1                                                           
         jm1    = j-1                                                           
      do 10161 i= lcyc1, lcyc2                                                  
         im1    = i-1                                                           
         ip1    = i+1                                                           
         if (im1.lt.1) im1 = l1                                                 
         if (ip1.gt.l) ip1 = 2                                                  
         r      = (x(ip1,j) + x(im1,j) - 2.*x(i,j))*zinv(j)                     
     &            + (x(i,jp1) + x(i,jm1) - 2.*x(i,j))*zzinv                     
         r      = (r-y(i,j))*z(j)                                               
         if (lsc-nsc) 29,29,30                                                  
   29    x(i,j) = x(i,j) + alfa*r                                               
   30    if (abs(r).le.eps)  nrel = nrel+1                                      
10161 continue                                                                  
      nsc       = nsc+1                                                         
      if (nrel-npts) 13,14,14                                                   
   14 if (lsc .ge. nsc) go to 300                                               
   18 lsc       = nsc+1                                                         
   13 if (nsc.lt.ia) go to 15                                                   
  300 continue                                                                  
      write (6,1000)                                                            
      write (6,1001)npts,nrel,nsc,ia                                            
 1000 format(20x,'progress of relaxation ntps nrel nsc ia')                     
 1001 format(6x,4i9)                                                            
      return                                                                    
      end                                                                       
c                                                                               
      subroutine STREAMF (u,v,ld,l  ,ll,m,dy,dx,a   ,psi)                            
c  this subroutine  computes the streamfunction field,                          
c  psi, from  the  equation, delsquared psi = relative                          
c  vorticity , given the horizontal components of  the                          
c  wind field.the sequential over-relaxation method is                          
c  used   . the  boundary conditions  imposed  for the                          
c  numerical  solution  of  the  equation  consist  of                          
c  specifying no net mass flux out of the domain.                               
c                                                                               
c  definitions    :                                                             
c                                                                               
c  u(ld,m)        : zonal wind component field.                                 
c  v(ld,m)        : meridional wind component field.                            
c  ld             : first  dimension  of arrays  u, v, a and                    
c                   psi.  if u, v  input  data is  cyclic in                    
c                   the zonal direction, set ld=l, if  it is                    
c                   not cyclic, set ld=ll.                                      
c  l              : number  of  grid  points  in  the  zonal                    
c                   direction   including  the  grid  points                    
c                   added to make domain cyclic in the zonal                    
c                   direction.                                                  
c  ll             : number  of  grid  points  in  the  zonal                    
c                   direction  excluding   the  grid  points                    
c                   added to make domain cyclic .                               
c  m              : number  of  grid  points  added  in  the                    
c                   meridional direction.                                       
c  dy             : grid size  in  meters  in the meridional                    
c                   direction.                                                  
c  dx(m)          : grid   size   in  meters  in  the  zonal                    
c                   direction for each row of grid points.                      
c  a(ld,m)        : relative vorticity field.                                   
c  psi(ld,m)      : streamfunction field.                                       
c                                                                               
      real u(ld,m),v(ld,m),psi(ld,m),a(ld,m),dx(m)                              
      real dxsq(100), dxsqinv(100)                                              
      data nsbd/999/                                                            
      l1        = l-1                                                           
      m1        = m-1                                                           
      l2        = l-2                                                           
      m2        = m-2                                                           
      ll1       = ll-1                                                          
      lcyc      = l-ll                                                          
      dysq      = dy*dy                                                         
      dysqinv   = 1./dysq                                                       
      do 10170 j= 1, m                                                          
      dxsq(j)   = dx(j)*dx(j)                                                   
      dxsqinv(j)= 1./dxsq(j)                                                    
10170 continue                                                                  
c                                                                               
c  compute the relative vorticity from the u and v field.                       
c                                                                               
      do 10171 i= 1, ld                                                         
      do 10171 j= 1, m                                                          
      psi(i,j)  = 0.                                                            
10171 continue                                                                  
      do 10173 j= 2, m1                                                         
      do 10172 i= 2, ld-1                                                       
10172 a(i,j)    = (v(i+1,j) - v(i-1,j))/(2.*dx(j))                              
     &              -(u(i,j+1) - u(i,j-1))/(2.*dy)                              
      a(1,j)    = (v(2,j) - v(ld-1,j))/(2.*dx(j))                                 
     &             -(u(1,j+1)-u(1,j-1))/(2.*dy)                                 
      a(ll,j)   = a(1,j)                                                        
10173 continue                                                                  
c                                                                               
c  adjusted the normal  outward component of the wind field                     
c  to   yield zero  outward  net mass flux.                                     
c                                                                               
      if (lcyc.eq.0) then                                                       
      vno1      = (v(1,m)+v(l,m))*dx(m)/2.                                      
      uno1      = (abs(v(1,m))+abs(v(l,m)))*dx(m)/2.                            
      do 10174 i= 2, l1                                                         
      uno1      = uno1 + abs(v(i,m))*dx(m)                                      
      vno1      = vno1 + v(i,m)*dx(m)                                           
10174 continue                                                                  
      vno2      = -(v(1,1) + v(l,1))*dx(1)/2.                                   
      uno2      = (abs(v(1,1)) + abs(v(l,1)))*dx(1)/2.                          
      do 10175 i= 2, l1                                                         
      vno2      = vno2 - v(i,1)*dx(1)                                           
      uno2      = uno2 + abs(v(i,1))*dx(1)                                      
10175 continue                                                                  
      eps1      = vno1/uno1                                                     
      eps2      = vno2/uno2                                                     
      do 10176 i= 1, l                                                          
      v(i,1)    = v(i,1) + eps2*abs(v(i,1))                                     
      v(i,m)    = v(i,m) - eps1*abs(v(i,m))                                     
10176 continue                                                                  
      else                                                                      
      vno       = v(1,m)*dx(m)/2. + v(ll,m)*dx(m)/2.                            
      uno       = abs(v(1,m))*dx(m)/2. + abs(v(ll,m))*dx(m)/2.                  
      do 10177 i= 2, ll1                                                        
      uno       = uno + abs (v(i,m))*dx(m)                                      
10177 vno       = vno + v(i,m)*dx(m)                                            
      vno       = vno + u(ll,m)*dy/2. + u(ll,1)*dy/2.                           
      uno       = uno + abs(u(ll,m))*dy/2. + abs(u(ll,1))*dy/2.                 
      do 10178 j= 2, m1                                                         
      uno       = uno + abs(u(ll,j))*dy                                         
10178 vno       = vno + u(ll,j)*dy                                              
      vno       = vno - v(ll,1)*dx(1)/2. - v(1,1)*dx(1)/2.                      
      uno       = uno + abs(v(ll,1))*dx(1)/2. + abs(v(1,1))*dx(1)/2.            
      do 10179 i= 2, ll1                                                        
      uno       = uno + abs(v(i,1))*dx(1)                                       
10179 vno       = vno - v(i,1)*dx(1)                                            
      vno       = vno - u(1,1)*dy/2. - u(1,m)*dy/2.                             
      uno       = uno + abs(u(1,1))*dy/2. + abs(u(1,m))*dy/2.                   
      do 10180 j= 2, m1                                                         
      uno       = uno + abs(u(1,j))*dy                                          
10180 vno       = vno - u(1,j)*dy                                               
      eps       = vno/uno                                                       
      do 10181 i= 1, ll                                                         
      v(i,1)    = v(i,1) + eps*abs(v(i,1))                                      
10181 v(i,m)    = v(i,m) - eps*abs(v(i,m))                                      
      do 10182 j= 1, m                                                          
      u(1,j)    = u(1,j) + eps*abs(u(1,j))                                      
10182 u(ll,j)   = u(ll,j) - eps*abs(u(ll,j))                                    
      endif                                                                     
c                                                                               
c  perform the sequential relaxation                                            
c                                                                               
      psi(1,m)  = 0.                                                            
      do 10183 i= 2, ll                                                         
10183 psi(i,m)  = psi(i-1,m) + (v(i,m) + v(i-1,m))*dx(m)/2.                     
      do 10184 jj= 1, m1                                                        
      j         = m-jj                                                          
10184 psi(ll,j) = psi(ll,j+1) + (u(ll,j) + u(ll,j+1))*dy/2.                     
      do 10185 ii= 1, ld                                                        
      i         = ld-ii + 1                                                        
      psi(i,1)  = psi(i+1,1) - (v(i,1) + v(i+1,1))*dx(1)/2.
10185 continue                     
      do 10186 j= 2, m1                                                         
10186 psi(1,j)  = psi(1,j-1) - (u(1,j) + u(1,j-1))*dy/2.                        
      if (lcyc.eq.0) then                                                       
      lx        = l                                                             
      x1        = l1                                                            
      else                                                                      
      lx        = ll                                                            
      lx1       = ll1                                                           
      endif                                                                     
      call relaxmod (psi,dysqinv,dxsq,a,dxsqinv,lx,lx1,m1,m,1)                  
      return                                                                    
      end                                                                       
c                                                                               
      subroutine zfieldmod (psi,upsi,vpsi,l, m, dx, dy,                         
     &                      cor,a2,z,fk1,fk2,zbar,beta)                         
c                                                                               
c  this subroutine  computes  the height field from the                         
c  streamfunction field by solving  the reverse balance                         
c  equation.  three scaled   version   of  the  reverse                         
c  balance equation are possible.  these are controlled                         
c  by the argument variables, fk1 and fk2. for,                                 
c  fk1 = 0  fk2 = 0 we have geostrophic zfield,                                 
c  fk1 = 1  fk2 = 0 we have linear balance zfield,                              
c  fk1 = 1  fk2 = 1 we have non-linear balance zfield.                          
c  other combinations are not allowed.                                          
c    note that this subroutine is similar in principles                         
c  to subroutine zfield where the cyclicity condition has                       
c  been added .                                                                 
c                                                                               
c  definitions   :                                                              
c                                                                               
c  psi(l,m)      : streamfunction field.                                        
c  upsi(l,m)     : rotational u component field.                                
c  vpsi(l,m)     : rotational v component field.                                
c  l             : number  of  grid  points  in  the  zonal                     
c                  direction.                                                   
c  m             : number of grid points  in the meridional                     
c                  direction.                                                   
c  dx(m)         : grid size in meters for each row of grid                     
c                  points in the zonal direction.                               
c  dy            : grid size  in  meters  in the meridional                     
c                  direction.                                                   
c  cor(m)        : coriolis parameter.                                          
c  a2(l,m)       : contains  the forcing function  field of                     
c                  the balance equation .                                       
c  z(l,m)        : height field.                                                
c  fk1,fk2       : control integers.  function as explained                     
c                  above.                                                       
c  zbar          : mean free surface height.                                    
c  beta          : beta plane parameter (=df/dy).                               
c                                                                               
      real  z(l,m),psi(l,m),dx(m),zinv(100),  z1(100)                           
      real  vpsi(l,m),  upsi(l,m),  a2(l,m),  cor(m)                            
c                                                                               
      l1        = l-1                                                           
      l2        = l-2                                                           
      m1        = m-1                                                           
      m2        = m-2                                                           
      ladd      = 6                                                             
      ginv      = 1.0/9.81                                                      
      do 10190 j= 1, m                                                          
      z1(j)     = dx(j)**2                                                      
10190 zinv(j)   = 1./z1(j)                                                      
      zz        = dy**2                                                         
      zzinv     = 1./zz                                                         
      sum       = 0.                                                            
      do 10191 i= 1, l                                                          
      do 10191 j= 1, m                                                          
      sum       = sum+psi(i,j)                                                  
10191 continue                                                                  
      psibar    = sum/(l*m)                                                     
c                                                                               
      do 10192 j= 2, m1                                                         
      do 10192 i= 1, l                                                          
      im1       = i-1                                                           
      ip1       = i+1                                                           
      if (im1.lt.1) im1 = l1                                                    
      if (ip1.gt.l) ip1 = 2                                                     
      vpsi(i,j) =  (psi(ip1,j) - psi(im1,j))/(2.0*dx(j))                        
      upsi(i,j) = -(psi(i,j+1)*dx(j)/dx(j+1) -                                  
     &            psi(i,j-1)*dx(j)/dx(j-1))/(2.*dy)                             
10192    continue                                                               
      call lapmod (z ,psi,dx,dy,l,m,1,l1,m1,l2,m2,ladd)                         
      call jacmod (a2,upsi,vpsi,dx,dy,l,m,1,l1,m1,l2,m2,ladd)                   
      do 10193 i= 1, l                                                          
      do 10193 j= 2, m1                                                         
      betau     = upsi(i,j)*beta                                                
      fvort     = cor(j)*z (i,j)                                                
      a2(i,j)   = fvort*ginv-betau*ginv*fk1+2.*a2(i,j)*ginv*fk2                 
10193 continue                                                                  
      do 10194 i= 1, l                                                          
      a2(i,1)   = 0.                                                            
      a2(i,m)   = 0.                                                            
10194 continue                                                                  
      do 10195 j= 1, m                                                          
      do 10195 i= 1, l1                                                         
      z(i,j)    = zbar + cor(j)*ginv*(psi(i,j) - psibar)                        
10195 continue                                                                  
      do 10196 j= 1, m                                                          
10196 z(l,j)    = z(1,j)                                                        
      call relaxmod (z,zzinv,z1,a2 ,zinv,l,l1,m1,m,0)                           
      return                                                                    
      end                                                                       
c                                                                               
      subroutine terr (h,l,m,del,slat,dhdy,dhdx,dx,hd)                          
c                                                                               
c  this subroutine computes zonal and meridional gradients                      
c  of the terrain fields .it outputs terrain (h),dhdx and                       
c  dhdy to tape23 for input to the shllow water model silepe .                  
c                                                                               
      real h(l,m),dhdy(l,m),dhdx(l,m),dx(m),hd(l,m)                             
      pi        = 4.*atan(1.)                                                   
      rad       = pi/180.                                                       
      dy        = del*111.1*1000.                                               
      do 10200 i= 1, m                                                          
      x         = (slat+del*(i-1))*rad                                          
      dx(i)     = dy*cos(x)                                                     
10200 continue                                                                  
      do 10201 j= 1, m                                                          
      do 10201 i= 1, l                                                          
      ip1       = i+1                                                           
      im1       = i-1                                                           
      if (i.eq.l) ip1 = 1                                                       
      if (i.eq.1) im1 = l                                                       
      dhdx(i,j) = (h(ip1,j)-h(im1,j))/(2.*dx(j))                                
10201 continue                                                                  
      do 10202 j= 1, m                                                          
      jp1       = j+1                                                           
      jm1       = j-1                                                           
      if (j.eq.1) jm1 = 1                                                       
      if (j.eq.m) jp1 = m                                                       
      dyy       = 2.0*dy                                                        
      if (j.eq.1) dyy = dy                                                      
      if (j.eq.m) dyy = dy                                                      
      do 10202 i= 1, l                                                          
      dhdy(i,j) = (h(i,jp1)-h(i,jm1))/(dyy)                                     
10202 continue                                                                  
      rewind (23)                                                               
      write(23,1000) ((h(i,j),i=1,l),j=1,m)                                     
      write(23,1000) ((dhdx(i,j),i=1,l),j=1,m)                                  
      write(23,1000) ((dhdy(i,j),i=1,l),j=1,m)                                  
 1000 format(6e13.6)                                                            
      return                                                                    
      end                                                                       
