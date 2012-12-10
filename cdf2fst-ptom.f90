
!
! NetCDF to Standard file convertor
! It can work in conjunction with the grib2cdf converter
!
! Andrew Ryzhkov and Michel Bourqui
! Atmospheric and Oceanic Sciences Department, McGill University
!

! ToDo:
!	* conversion units of the fields: GRIB/NetCDF => FST
!	* variable parameters store in a file
!	* level parameters store in a file
!	* pressure and sigma levels

program cdf2fst

  implicit none

  include 'netcdf.inc'
   
  type FSTVar
     character*64 	longname
     character*32 	cdfname
     character*4  	fstname
     character*8 	units
     real		factor
     real		bias
     integer		cdfid
  end type FSTVar

!  integer, parameter :: VarCnt = 6
!  type (FSTVar) :: Vars(VarCnt)
!  data Vars /	&
!	FSTVar( "Surface pressure ","PS","P0", "mbar ", 1.00000,   0.00, 0 ),	&
!	FSTVar( "Geopotential     ","Z ","GZ", "dm   ", 0.10000,   0.00, 0 ),	&	! convert from meters to decameters *0.1
!	FSTVar( "Air Temperature  ","T ","TT", "deg C", 1.00000,   0.00, 0 ),	&	! for Kelvin use 273.15
!	FSTVar( "Zonal wind       ","U ","UU", "knots", 1.94384,   0.00, 0 ),	&	! convert from m/2
!	FSTVar( "Meridional wind  ","V ","VV", "knots", 1.94384,   0.00, 0 ),	&	! to knots
!	FSTVar( "Specific humidity","Q ","HU", "kg/kg", 1.00000,   0.00, 0 )	/

!! MOZART4 variables
!  integer, parameter :: VarCnt = 8
!  type (FSTVar) :: Vars(VarCnt)
!  data Vars /	&
!	FSTVar( "Surface pressure"        ,"PS"           ,"P0"   , "Pa"   , 0.01000 ,   0.00, 0 ) ,   &               ! Pa => mbar
!	FSTVar( "Hydroxyl radical"        ,"OH_VMR_inst"  ,"OH"   , "kg/kg", 0.587148,   0.00, 0 ) ,   &               ! VMR*Moh/Mair => MMR
!        FSTVar( "Peroxyl radical"         ,"HO2_VMR_inst" ,"HO2"  , "kg/kg", 1.139499,   0.00, 0 ) ,   &
!        FSTVar( "Ozone"                   ,"O3_VMR_inst"  ,"O3"   , "kg/kg", 1.657053,   0.00, 0 ) ,   &  
!        FSTVar( "Hydrogen peroxide"       ,"H2O2_VMR_inst","H2O2" , "kg/kg", 1.174297,   0.00, 0 ) ,   &  
!        FSTVar( "Nitrogen dioxide"        ,"NO2_VMR_inst" ,"NO2"  , "kg/kg", 1.588259,   0.00, 0 ) ,   &  
!        FSTVar( "Sulfur dioxide"          ,"SO2_VMR_inst" ,"SO2"  , "kg/kg", 2.211766,   0.00, 0 ) ,   &  
!        FSTVar( "Black carbon hydrofilic" ,"CB2_VMR_inst" ,"SOOT" , "kg/kg", 0.414648,   0.00, 0 )        /  

! p-TOMCAT variables
  integer, parameter :: VarCnt = 4
  type (FSTVar) :: Vars(VarCnt)
  data Vars /   &
!        FSTVar( "Pressure"        ,"p"   ,"P"   , "mbar" , 0.010000    ,   0.00, 0 ) ,   &               ! Pa => mbar
        FSTVar( "Nitrogen dioxide","NO2" ,"NO2" , "kg/kg", 1.588259e-9 ,   0.00, 0 ) ,   &               ! ppbv => mmr
        FSTVar( "Hydroxyl radical","OH"  ,"OH"  , "kg/kg", 0.587148e-12,   0.00, 0 ) ,   &               ! pptv*Moh/Mair => mmr
        FSTVar( "Bromine"         ,"Br"  ,"BR"  , "kg/kg", 2.758406e-12,   0.00, 0 ) ,   &
        FSTVar( "Bromine oxide"   ,"BrO" ,"BRO" , "kg/kg", 3.310778e-12,   0.00, 0 )        /

  ! FST variables
  character*4 nomvar
  character*1 typvar, grtyp
  character*8 etiket
  
  integer :: dateo, datev, deet, npas, ni, nj, nk, npak, datyp, yyyymmdd, gdin
  integer :: ip1, ip2, ip3, ip1z, ip2z, ip3z
  integer :: ig1, ig2, ig3, ig4
  integer, dimension(2) :: fullDate

  ! RPN FST functions externals
  integer, external :: fnom, fstouv, fclos, fstfrm, newdate, fstecr, write_encode_hyb, hybref_to_ig, read_decode_hyb, ezqkdef, gdll

  ! other variables
  integer ier,i,iun,xi,yi,ii,jj,kk
  real work(120, 32)
  real, allocatable, dimension (:) :: aklay, bklay

  ! NetCDF variables
  character*256, allocatable, dimension(:)     :: dimnames, varnames			! dimensions and variables names
  integer,       allocatable, dimension(:)     :: dimlen, vartype, vardim, varnatt	! dimensions lengths, variables types, variables dimensions, variables attributes
  integer,       allocatable, dimension(:,:)   :: vardimids
  real*8 ,       allocatable, dimension(:)     :: timearr				! times
  real   ,       allocatable, dimension(:)     :: lon, lat                              ! lon, lat arrays
  real   ,       allocatable, dimension(:,:,:) :: datarr, fstarr			! data
  real   ,       allocatable, dimension(:,:)   :: ps, zlat, zlon					! surface pressure
  character*256 inpfile, outfile, arg, title, datetimestr, datetimestr0
   
  integer stat, ndim, nvar, natt, unlimdid, nicid, ivar
  integer nlev, start(4), count(4)
  integer ilon, ilat, ilev, itim
  integer lons, lats, levs, tims, lev
  integer time

  ! hybrid levels parameters
  logical            :: hybrid = .true.							! .false.
  integer, parameter :: nkhyb = 28 !60                                                  ! levels
  real   , parameter :: Pref  = 800.							! reference pressure
  real   , parameter :: Ptop  = 10.0							! top pressure
  real   , parameter :: Rcoef = 1.0							! R coefficient
  real, dimension(nkhyb) :: hyb, phyb

  real, dimension(60) :: hyb60
  data hyb60 /								&		! 60 hybrid levels
       1.000000, 0.995992, 0.987042, 0.973284, 0.954921, 0.932221,	&
       0.905510, 0.875166, 0.841608, 0.805288, 0.766678, 0.726264,	&
       0.684534, 0.641967, 0.599026, 0.557492, 0.517476, 0.479070,	&
       0.442346, 0.407360, 0.374149, 0.342735, 0.313125, 0.285309,	&
       0.259269, 0.235802, 0.214639, 0.195536, 0.178280, 0.162679,	&
       0.148564, 0.135783, 0.124200, 0.113695, 0.104160, 0.094973,	&
       0.086184, 0.077831, 0.069945, 0.062546, 0.055649, 0.049257,	&
       0.043369, 0.037977, 0.033068, 0.028623, 0.024622, 0.020949,	&
       0.017617, 0.014628, 0.011976, 0.009648, 0.007626, 0.005889,	&
       0.004410, 0.003165, 0.002127, 0.001269, 0.000568, 0.000000 	&
  /

  !real, dimension(28) :: hyb28
  data hyb /						&				! 28 hybrid levels
    1.000, 0.993, 0.980, 0.955, 0.922, 0.884, 0.842, 0.796, &
    0.744, 0.688, 0.631, 0.574, 0.516, 0.460, 0.405, 0.351, &
    0.302, 0.258, 0.219, 0.185, 0.155, 0.127, 0.101, 0.075, &
    0.051, 0.027, 0.011, 0.000  /

  ! pressure level parameters
!  real, dimension(nkhyb) :: pres
!  data pres /						&
!    1000, 993, 988, 981, 975, 962, 950,                 &
!     937, 925, 913, 900, 875, 850, 825, 		&
!     800, 775, 750, 700, 650, 600, 550,			&
!     500, 400, 300, 200, 100,  50,  10  /

  ! and variables
  integer :: k
  real, dimension(nkhyb) :: A,B
  !real :: Pr

! from the GEM sources - the hybrid_to_pres function:
! 	hybm_8(k)= hyb(k) + (1-hyb(k)) * ptop/pref
! 	prpref = 100.*ptop/hybm_8(1)
! 	pr1 = 1./(1. - hybm_8(1))				! hybm_8(1) = ptop/pref
! 	pibb(k)  = ((hybm_8(k) - hybm_8(1))*pr1 ) ** rcoef
! 	pia(k)  = prpref * ( hybm_8(k) - pibb(k) )
! 	pressure(i,k) = (pia(k)+pibb(k)*ps(i))
! follows that A and B are
! 	B(k) = hyb(k)**rcoef, 
!	A(k) = Pref*(hyb(k)-B(k)) + Ptop*(1-hyb(k))
!	P(i) = A(k) + B(k)*Ps(i)
  
  ! arguments variables
  integer iargc
  
  
  if (hybrid) then
    do k=1,nkhyb									! calculate A and B parameters for the hybrid levels
      B(k) = hyb(k)**Rcoef
      A(k) = Pref*(hyb(k)-B(k)) + Ptop*(1-hyb(k))
      !print '(3f12.6)',hyb(k),A(k),B(k)
    end do
  endif
!  print '(6(f9.6,","))', (hyb(k),k=nkhyb,1,-1)

  ! initial file names
  inpfile = ''
  outfile = ''

  ! Arguments
  if ( iargc()==0 ) call usage_stop
  
  i=1
  do while (iargc().ge.i)
       call getarg(i,arg)
       i=i+1
       select case (arg)
       case ('-i')
           call getarg(i,inpfile)
           i=i+1
       case ('-o')
           call getarg(i,outfile)
           i=i+1
       case ('-h')
           call usage_stop
       end select
  enddo
   
  if ( len_trim(inpfile) == 0 ) then
    print *, 'No input file specified with "-i" argument'
    call usage_stop
  end if

  if ( len_trim(outfile) == 0 ) then
    print *, 'No output file specified with "-o" argument'
    call usage_stop
  end if

  ! Program start
  print '("Open input file: <",a,"> ",$)', trim(inpfile)
  stat = nf_open( trim(inpfile), nf_nowrite, nicid )
  call handle_err( stat, "Open input file" )
  print *,'- Ok.'
      
  ! Read header
  print '(a,$)',"Read header "
  stat = nf_inq( nicid, ndim, nvar, natt, unlimdid )
  call handle_err( stat, "Read header" )
  print *,'- Ok.'
   
  !stat = nf_inq_attlen  ( nicid, NF_GLOBAL, 'title', i )
  !stat = nf_get_att_text( nicid, NF_GLOBAL, 'title', title )
  !call handle_err( stat, "Get title attribute" )
  !print *,'Title: <',TRIM(title(1:i)),'>'

  print '(a,i8)',     &
     'Dimensions:',ndim, &
     'Variables :',nvar, &
     'Attributes:',natt, &
     'Unlimited :',unlimdid
   
  ! Read dimensions information
  print '(a,$)','Read dimensions '
  allocate( dimnames(ndim) )
  allocate( dimlen  (ndim) )
  do i=1,ndim
     stat=nf_inq_dim(nicid,i,dimnames(i),dimlen(i))
     call handle_err( stat, "Read dimensions" )
  end do
  print *,'- Ok.'

  print '("(",i2,")",a12,":",i8)',((i,trim(dimnames(i)),dimlen(i)),i=1,ndim)
   
  do i=1,ndim
    if ( trim(dimnames(i)) == 'lon') then
      ilon = i
      lons = dimlen(i)
    endif
    if ( trim(dimnames(i)) == 'lat') then
      ilat = i
      lats = dimlen(i)
    endif
    if ( trim(dimnames(i)) == 'niv') then
      ilev = i
      levs = dimlen(i)
    endif
    if ( trim(dimnames(i)) == 'time') then
      itim = i
      tims = dimlen(i)
    endif
  enddo
  
  ! Read variables information
  print '(a,$)','Read variables '
  allocate(varnames(nvar))
  allocate(vartype(nvar))
  allocate(vardim(nvar))
  allocate(vardimids(nvar,ndim))
  allocate(varnatt(nvar))
  allocate(lon(lons+1),lat(lats))
  do i=1,nvar
     stat=nf_inq_var(nicid,i,varnames(i),vartype(i),      &
                     vardim(i),vardimids(i,:),varnatt(i))
     call handle_err(stat,"Read variables")
  end do
  print *,'- Ok.'
   
  allocate(timearr(tims))
  allocate(aklay(levs),bklay(levs))

  !print '(a10,3i8)',((varnames(i),vartype(i),vardim(i),varnatt(i)),i=1,ndim)
  !print '(7i6)', ((vardimids(i,j),j=1,ndim),i=1,nvar)
  do ivar=1, nvar
  
    print '("(",i3,")",a20,"[1:",i1,"] depend on: ",8a8)', ivar,trim(varnames(ivar)), vardim(ivar), (trim(dimnames(vardimids(ivar,i))),i=1,vardim(ivar))
    
    do i=1,VarCnt
      if ( TRIM(Vars(i)%cdfname) == TRIM(varnames(ivar)) ) Vars(i)%cdfid=ivar	! assign NetCDF variables IDs
    enddo
    
    if ( 'vca' == TRIM(varnames(ivar)) ) then
      stat=nf_get_var_real ( nicid, ivar, aklay ) 
      call handle_err(stat,"Read aklay variable")
    endif
    
    if ( 'vcb' == TRIM(varnames(ivar)) ) then
      stat=nf_get_var_real ( nicid, ivar, bklay ) 
      call handle_err(stat,"Read bklay variable")
    endif
    
    if ( 'lon' == TRIM(varnames(ivar)) ) then
      stat=nf_get_var_real ( nicid, ivar, lon ) 
      call handle_err(stat,"Read lon variable")
    endif
    
    if ( 'lat' == TRIM(varnames(ivar)) ) then
      stat=nf_get_var_real ( nicid, ivar, lat ) 
      call handle_err(stat,"Read lat variable")
    endif
    
    if ( 'time' == TRIM(varnames(ivar)) ) then
      stat = nf_inq_attlen  ( nicid, ivar, 'units', i )
      call handle_err( stat, "Get time attribute" )
      stat = nf_get_att_text( nicid, ivar, 'units', datetimestr )
      call handle_err( stat, "Get time attribute" )
      ! read initial time and convert them to FST format
      print *,'Long   time: <',datetimestr(1:i),'>'
      write (datetimestr0,'(a4,2a2)') datetimestr(12:15),datetimestr(17:18),datetimestr(20:21)
      print *,'Short  time: <',trim(datetimestr0),'>'
      read  (datetimestr0,'(i)') yyyymmdd
      print '(a,i8,a)','Integer time: <',yyyymmdd,'>'
      stat = nf_get_var_double( nicid, ivar, timearr )
      print '("Times: ",320f10.2)',timearr
      call handle_err(stat,"Read time variable")
    endif
    
    !if ( 'P0' == TRIM(varnames(ivar)) ) then 
    !  stat=nf_get_var_real ( nicid, ivar, Pr ) 
    !  call handle_err(stat,"Read P0 variable")
    !endif
    
  end do

  lat(1:lats) = lat(lats:1:-1)                                         ! reverse lat from 90..-90 to -90..90
  aklay = aklay(levs:1:-1)*1000.                         ! kPa => mbar
  bklay = 0     !bklay(levs:1:-1)
  do lev=levs-1,1,-1
    aklay(lev+1) = ( aklay(lev) + aklay(lev+1) )/2.       ! interface to middle
    !bklay(lev) = ( bklay(lev) + bklay(lev+1) )/2.
  enddo
!  do lev=levs,1,-1
!    bklay(lev+1) = bklay(lev)
!  end do
!  bklay(1) = 1.0
!  print '("HYAM=",32f10.4)',aklay(1:levs)
!  print '("HYBM=",32f10.4)',bklay(1:levs)
  print '("lon(",i3,")=",320f7.2)',lons,lon
  print '("lat(",i3,")=",320f7.2)',lats,lat

!  print '(a,i3,a,a,a,i3,a)', "Variable#", ivar, " '", trim(varnames(ivar)), "', which has ", vardim(ivar), " dimensions"
!  print '(a,4i12)', "Namely   : ", (vardimids(ivar,i),i=1,4) 
!  print '(a,4i12)', "with dim : ", (dimlen(vardimids(ivar,i)),i=1,4)
!  print '(a,4a12)', "and names: ", (trim(dimnames(vardimids(ivar,i))),i=1,4)
   
!  print *, "Allocate ",dimlen(vardimids(ivar,1))*dimlen(vardimids(ivar,2))*4," bytes"
  allocate( datarr( lons  , lats, levs  ), stat=stat );	if (stat /= 0) stop 'Out of memory'
  allocate( fstarr( lons+1, lats, nkhyb ), stat=stat );	if (stat /= 0) stop 'Out of memory'
!  allocate( gz    ( lons+1, lats        ), stat=stat );	if (stat /= 0) stop 'Out of memory'
  allocate( ps    ( lons+1, lats        ), stat=stat );	if (stat /= 0) stop 'Out of memory'
   
!  allocate(fstarr(dimlen(vardimids(ivar,1))+1, dimlen(vardimids(ivar,2))), stat=stat )		! 361 x 181
!  if (stat /= 0) stop 'Out of memory'

  ! Association of the RPN standard file produced by the program with the FORTRAN logical unit 1.
  iun = 1
  ier = fnom(iun, outfile, 'STD+RND', 0)
  if (ier<0) then
     print *, 'Fatal error while opening the file (FNOM)'
     stop
  endif

  ! Open the FST file
  iun = 1
  ier = fstouv(iun, 'RND')
  if ( ier<0 ) then
     print *, 'Cannot open unit:', iun,' in random access mode (FSTOUV)'
     stop
  endif

  ip2 = 0
  ip3 = 0

  ! Initialization of the standard file attributes that remain constant for all fields
  typvar = 'P'
  etiket = 'p-TOMCAT'

  ip1 = 0
  ip2 = 0
  ip3 = 0

  ni = lons+1;          lon(lons+1)=lon(lons)+(lon(2)-lon(1))
  nj = lats
  nk = 1

  deet = 0
  npas = 0
  dateo = 0
  
  datyp =  1
  npak  = -32   !-24   !-32   !-16

  grtyp  = 'L'								! grid type:	'L' - cylindrical equidistant (lat-lon).
  !call cxgaig(grtyp, ig1, ig2, ig3, ig4, -89.5, -179.5, 1., 1.)		! convert grid parameters
  !lat(1), lon(1), lat(2)-lat(1), lon(2)-lon(1)
  !call cxgaig(grtyp, ig1, ig2, ig3, ig4, lat(1), lon(1), lat(2)-lat(1), lon(2)-lon(1) )         ! convert grid parameters
  call cxgaig(grtyp, ig1, ig2, ig3, ig4, 0., 0., 1., 1. )         ! convert grid parameters
                                                                                        !XLAT0: latitude of the southwest corner of the grid.
                                                                                        !XLON0: longitude of the southwest corner of the grid.
                                                                                        !DLAT: latitudinal grid length in degrees.
                                                                                        !DLON: longitudinal grid length in degrees.
  ip1z = 4200
  ip2z = 4201
  ip3z = 4202
  
  nomvar='>>'
  ier = fstecr(lon, lon,             &
               npak, iun, dateo, deet, npas,             &
               ni, 1, 1,                                 &
               ip1z, ip2z, ip3z, typvar, nomvar, etiket, grtyp,       &
               ig1, ig2, ig3, ig4, datyp, .false. )
  
  nomvar='^^'
  ier = fstecr(lat, lat,             &
               npak, iun, dateo, deet, npas,             &
               1, nj, 1,                                 &
               ip1z, ip2z, ip3z, typvar, nomvar, etiket, grtyp,       &
               ig1, ig2, ig3, ig4, datyp, .false. )
               
  grtyp  = 'Z'                                                          ! grid type:    'Z' - cartesian grid with a non-constant mesh
  gdin    = ezqkdef(ni, nj, grtyp, ip1z, ip2z, ip3z, 0, iun)
  
  allocate(zlat(ni,nj))
  allocate(zlon(ni,nj))
  ier     = gdll(gdin, zlat, zlon)
  print '("lon(",i3,")=",320f7.2)', lons, zlon(:,1)
  print '("lat(",i3,")=",320f7.2)', lats, zlat(1,:)
  nomvar = 'LA'
  ier=fstecr(zlat, zlat, npak, iun, dateo, deet, npas, ni, nj, nk, &
       ip1, ip2, ip3, typvar, nomvar, etiket, grtyp, ip1z, ip2z, ip3z, 0, datyp, .false.)
  nomvar = 'LO'
  ier=fstecr(zlon, zlon, npak, iun, dateo, deet, npas, ni, nj, nk, &
       ip1, ip2, ip3, typvar, nomvar, etiket, grtyp, ip1z, ip2z, ip3z, 0, datyp, .false.)
  deallocate(zlat,zlon)
  
  ig1=ip1z;  ig2=ip2z;  ig3=ip3z
  
  
  !yyyymmdd = 20050101                                                           ! 732312 days = 20050101
  !timearr  = timearr  - 732312.0                                                ! change units from "days since 0-0-0" to "days since 2005-01-01"
  ier = newdate(dateo, yyyymmdd, 0, 3)                                          ! obtain date
  datev = dateo

  ! Main cycle over time periods
  do time=1, 2!tims							! loop over time

    if ( tims == 1 ) then
      deet = 6*3600
    else
      deet = (timearr(2)-timearr(1))*24*3600				! suggest constant time step, in sec
    endif
    npas = timearr(time)*24*3600/deet						! step number
    ip2 = ((npas*deet+1800)/3600)
    !print '("ip2,deet,npas:",3i10)',ip2,deet,npas
                                                                                                                           
    call incdat(datev, dateo, (deet*npas+1800)/3600 )                                                                      
    ier = newdate ( datev, fullDate(1), fullDate(2), -3 )
    
    print ('"Time: ",i6,f8.2," => ",i8.8," ",i8.8'), time, timearr(time), fullDate
    
    if ( time == 1 ) then
      if ( hybrid ) then                                                                                     !
        ier = write_encode_hyb (iun,'HY',ip2,ip3,etiket,datev,	&	! encode and write information about hybrid levels !
                            ptop,pref,rcoef)                                                                                 !
        !print *,"write_encode_hyb: ",ier                    
      else
        fstarr(1:lons+1, 1:lats, 1) = Ptop				! top pressure
        ier = fstecr(fstarr(1:lons+1, 1:lats, 1),		&
                    WORK, npak, iun, datev, 0, 1,		&
                    ni, nj, nk,					&
                    0, 0, 0, typvar, 'PT', etiket, grtyp,	&
                    ig1, ig2, ig3, ig4, datyp, .false. )
      endif
    endif                                                                                                  !
    
    do ivar=1, nvar
      if ( 'p' == TRIM(varnames(ivar)) ) exit                         ! find cdfid for pressure variable 'p'
    enddo
!    do i=1, VarCnt                                                      ! loop over variables to find pressure id
!        if ( trim(Vars(i)%cdfname) == 'p' ) then   
!            ivar = Vars(i)%cdfid
!            ii = i
!        endif
!    enddo
    
    start=(/    1,    1,    1, time /)  
    count=(/ lons, lats, levs, 1    /)
    stat=nf_get_vara_real ( nicid, ivar, start, count, datarr )
    call handle_err(stat,"Read pressure")
    fstarr(1:lons,1:lats,1:levs) = datarr(1:lons,1:lats,levs:1:-1)*0.01      ! revert verticaly and convert units
    fstarr(lons+1,1:lats,1:levs) = fstarr(1,1:lats,1:levs)                                              ! longitude copy: 0=>360
    
    ps(1:lons+1,1:lats) = fstarr(1:lons+1,1:lats,1)                                                     ! surface pressure
    !print *, ps
    !print *,fstarr(lons+1,1:lats,1)
    print *,"Write surface pressure"
    !print *,ps(lons+1,1:lats)
    call convip ( ip1, hyb(nlev), 1, 2, etiket, .false. )                                               ! pressure level
    ier = fstecr(ps(1:lons+1, lats:1:-1),                          &
                work, npak, iun, dateo, deet, npas,                &
                ni, nj, nk,                                        &
                ip1, ip2, ip3, typvar, 'P0', etiket, grtyp,        &
                ig1, ig2, ig3, ig4, datyp, .false. )
    
    ! Calculate B for mid-points from A and P
    do lev=levs,1,-1
        
        !print *,ps(1:2,1:2)
        !ps = datarr(1:lons,1:lats,1)
        !ps = ps*Vars(i)%factor + Vars(i)%bias                               ! convert units
        !print *, ps(lons+1,1:lats)
        !print *, ps(1,1:lats)
        !ps(lons+1,1:lats) = ps(1,1:lats)                                    ! longitude copy: 0=>360
        !print *, ps(lons+1,1:lats)
        bklay(lev) = sum   ( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ) / dble(lons*lats)
        
        !print '(i3,5f14.8)', lev, aklay(lev), &
        !                          sum   ( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ) / dble(lons*lats), &
        !                          minval( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ),&
        !                          maxval( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ),&
        !                          sum     ( dble(fstarr(1:lons,1:lats,lev))                                                ) / dble(lons*lats)
    enddo
    
    print '("HYAM=",32f10.4)',aklay(1:levs)
    print '("HYBM=",32f10.4)',bklay(1:levs)
    
!    do lev=1,levs
!      do xi=1,lons
!        do yi=1,lats
!          !fstarr(xi,yi,lev) = ( ( aklay(lev) + bklay(lev)*ps(xi,yi) ) + ( aklay(lev+1) + bklay(lev+1)*ps(xi,yi) ) )/2
!          fstarr(xi,yi,lev) = aklay(lev) + bklay(lev)*ps(xi,yi)
!        enddo
!      enddo
!      call convip ( ip1, hyb(lev), 1, 2, etiket, .false. )               ! level
!      fstarr(lons+1,1:lats,lev) = fstarr(1,1:lats,lev)
!      ier = fstecr( fstarr(1:lons+1, lats:1:-1, lev),                  &
!                    work, npak, iun, dateo, deet, npas,                &
!                    ni, nj, nk,                                        &
!                    ip1, ip2, ip3, typvar, 'PP', etiket, grtyp,        &
!                    ig1, ig2, ig3, ig4, datyp, .false. )
!    enddo

    do i=1, VarCnt							! loop over variables

      ivar = Vars(i)%cdfid
      
      !if ( trim(Vars(i)%cdfname) == 'PS' ) then
      !  ps   = datarr ( 1:lons, 1:lats, 1 )                     ! surface pressure => ps
      !  nlev = 1
      !else
        nlev = levs
      !endif
      
      start=(/    1,    1,    1,   time /)  
      count=(/ lons, lats, nlev,   1    /)

      print ('a16,"(",i3,")"'),Vars(i)%cdfname,nlev

      stat=nf_get_vara_real ( nicid, ivar, start, count, datarr )
      call handle_err(stat,"Read variable: "//Vars(i)%cdfname)

!      datarr = cshift ( datarr, -180 )					! shift coordinates to Greenwich (1:lons, 1:lats, 1:lev)
      if ( nlev > 1) datarr = datarr(:,:,nlev:1:-1)			! revert verticaly
      datarr = datarr*Vars(i)%factor + Vars(i)%bias			! convert units

!      if ( trim(Vars(i)%cdfname) == 'p' ) then
!        ps   = datarr ( 1:lons, 1:lats, 1 )			! surface pressure => ps
!        !nlev = 1
!      endif  

      !if ( .false. ) then
      if ( nlev/=1 ) then

        if ( hybrid ) then
!firstprivate(lons,lats,levs) 
!private(pp,phyb,xx,ier) 
!shared(datarr,aklay,bklay,ps,pr,hybrid,A,B,hyb,nkhyb,fstarr) 
!$omp parallel 
!$omp do 
          do xi=1,lons
            do yi=1,lats
              call CSpInt ( levs , aklay + bklay*ps(xi,yi)  , datarr(xi,yi,1:levs),           &
                            nkhyb, A + B*ps(xi,yi)          , fstarr(xi,yi,1:nkhyb)  )
            end do
          end do
!$omp end do
!$omp end parallel
        else
!$omp parallel 
!$omp do 
          do xi=1,lons
            do yi=1,lats
              call CSpInt ( levs , aklay + bklay*ps(xi,yi)     , datarr(xi,yi,1:levs),           &
                            nkhyb, Ptop*(1-hyb) + hyb*ps(xi,yi), fstarr(xi,yi,1:nkhyb)  )
            end do
          end do
!$omp end do
!$omp end parallel
        end if
        nlev = nkhyb
      
      else
        fstarr(1:lons,1:lats,1:nlev) = datarr(1:lons,1:lats,1:nlev)
      endif

      fstarr(lons+1,1:lats,1:nlev) = fstarr(1,1:lats,1:nlev)		! longitude copy: 0=>360
      
      do lev=nlev, 1, -1
        
!        if ( trim(Vars(i)%cdfname) == 'p' ) then
!            print '(i3,5f14.8)', lev, aklay(lev), &
!                                      sum   ( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ) / dble(lons*lats), &
!                                      minval( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ),&
!                                      maxval( ( dble(fstarr(1:lons,1:lats,lev)) - dble(aklay(lev)) ) / dble(ps(1:lons,1:lats)) ),&
!                                      sum     ( dble(fstarr(1:lons,1:lats,lev))                                                ) / dble(lons*lats)
!        endif

!        if ( trim(Vars(i)%cdfname) == 'p' .and. lev==1 ) then
!            call convip ( ip1, hyb(nlev), 1, 2, etiket, .false. )	! surface pressure
!            print *,"Write surface pressure"
!            ! Write a standard file record
!            ier = fstecr(fstarr(1:lons+1, lats:1:-1, lev:lev),                   &
!                        work, npak, iun, dateo, deet, npas,                &
!                        ni, nj, nk,                                        &
!                        ip1, ip2, ip3, typvar, 'P0', etiket, grtyp,        &
!                        ig1, ig2, ig3, ig4, datyp, .false. )
!            
!        else
            !print *,"Write "//Vars(i)%cdfname//" as "//Vars(i)%fstname//" at level ", lev
            !if (hybrid) then
              call convip ( ip1, hyb(lev), 1, 2, etiket, .false. )	! convert hybrid coordinate level to ip1
            !else
            !  call convip ( ip1, hyb(lev), 1, 2, etiket, .false. )	! convert sigma level to ip1
            !endif
!        endif
      
        nomvar = Vars(i)%fstname

!        do kk=1,nlev
!          do ii=1,lons+1
!            do jj=lats/2+1,1,-1
!              if ( fstarr(ii,jj,kk)<0. .or. fstarr(ii,jj,kk)>10000. ) then	! missing values
!                fstarr(ii,jj,kk) = fstarr(ii,jj+1,kk)
!              endif
!            enddo
!            do jj=lats/2,lats
!              if ( fstarr(ii,jj,kk)<0. .or. fstarr(ii,jj,kk)>10000. ) then	! missing values
!                fstarr(ii,jj,kk) = fstarr(ii,jj-1,kk)
!              endif
!            enddo
!          enddo
!        enddo

        where ( fstarr < 0. ) fstarr = 0.                                             ! remove negative values

        ! Write a standard file record
        ier = fstecr(fstarr(1:lons+1, lats:1:-1, lev:lev),	 		&
                     work, npak, iun, dateo, deet, npas,		&
                     ni, nj, nk,					&
                     ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,	&
                     ig1, ig2, ig3, ig4, datyp, .false. )
      enddo	! lev
    
    enddo	! i

  enddo		! time
  
  ! Close the standard file
  ier = fstfrm(1)

  ! Unlink the unit 1 from the file 
  ier = fclos(1)

  ! Close the file
  stat = nf_close( nicid )
  call handle_err( stat, "Close NetCDF file" )

  ! Deallocate memory
  deallocate (dimnames)
  deallocate (dimlen)
  deallocate (varnames)
  deallocate (vartype)
  deallocate (vardim)
  deallocate (vardimids)
  deallocate (varnatt)
  deallocate (fstarr)
  deallocate (datarr)
  deallocate (ps)
!  deallocate (gz)
  deallocate (timearr)
  deallocate(lat,lon)
  deallocate(aklay,bklay)
  
  stop 'Finish'
  

contains


  subroutine handle_err ( iret, msg )
    implicit none
    integer, intent(in) ::  iret
    character(len=*), optional, intent(in) :: msg
    include 'netcdf.inc'
    if ( iret /= NF_NOERR ) then
      if (present(msg)) then 
        print *, "NetCDF error: ", msg, nf_strerror(iret)
      else
        print *, "NetCDF error: ",      nf_strerror(iret)
      endif
      stop
    endif
  end subroutine handle_err
  

  pure integer function LinInt (xn,x,y,xin,xi,yi)
  ! Levels linear interpolation

  implicit none
  
  real, dimension(:), intent (in)       :: x,y,xi
  integer,            intent (in)       :: xn,xin
  real, dimension(:), intent (out)      :: yi
  real, dimension(xn-1)                 :: a,b
  real                                  :: ai, bi
  integer                               :: i,j

    !yi(1:xin) = y(1:xin); return

    do i=1, xn-1
      a(i) =   ( y(i)-y(i+1) ) 	           / ( x(i)-x(i+1) )
      b(i) = - ( x(i+1)*y(i)-x(i)*y(i+1) ) / ( x(i)-x(i+1) )
    end do
    
    !yi(1:xin) = y(1:xin); !return
    
    do i=1, xin

      if ( xi(i)>=x(1) ) then
!        ai = a(1)
!        bi = b(1)
        yi(i) = y(1)							! boundary value
      else
        if ( xi(i)<=x(xn) ) then
!          ai = a(xn-1)
!          bi = b(xn-1)
           yi(i) = y(xn)						! boundary value
        else
          do j=1, xn-1
            if ( xi(i)<=x(j) .and. xi(i)>=x(j+1) ) then
              ai = a(j)
              bi = b(j)
              exit
            end if
          end do
          yi(i) = ai*xi(i) + bi
        endif
      endif

    end do
    
    LinInt = 0

  end function LinInt
  

  pure subroutine CSpInt (xn,x,y,xin,xi,yi)
  ! Cubic-spline approximation

  implicit none
  
  real, dimension(xn) , intent (in)  :: x,y
  integer             , intent (in)  :: xn,xin
  real, dimension(xin), intent (in)  :: xi
  real, dimension(xin), intent (out) :: yi
  real, dimension((xn+1))            :: p2
  real :: xx,dx,a,b,c,d,x1,x2

  integer :: i,k

    call Cubic_Spline( xn-1, x(xn:1:-1), y(xn:1:-1), p2 )

    do i=1, xin

      xx = xi(xin-i+1)

      ! Find the interval that x resides
      if ( xx<=x(xn) ) then
        dx = xx-x(xn)
        k  = 1
      else
        if ( xx>=x(1) ) then
          dx = xx-x(1)
          k  = xn-1
        else
          k = 1
          dx = xx-x(xn)
          DO WHILE (dx>=0)
            k = k + 1
            dx = xx-x(xn-k+1)
          END DO
          k = k - 1
        endif
      endif

      ! Find the value of function f(x)
      dx = x(xn-(k+1)+1) - x(xn-K+1)
      a  =  p2(K+1)/(6*dx)
      b  = -p2(K  )/(6*dx)
      c  = -dx*p2(k+1)/6 + y(xn-(k+1)+1)/dx
      d  =  dx*p2(k  )/6 - y(xn- k   +1)/dx
      x1 = (Xx-X(xn-K    +1))
      x2 = (Xx-X(xn-(K+1)+1))
      yi(xin-i+1) = a*x1**3 + b*x2**3 + c*x1 + d*x2

    end do
    
  end subroutine CSpInt
  

pure SUBROUTINE CUBIC_SPLINE ( n, XI, FI, P2 )
  ! Function to carry out the cubic-spline approximation
  ! with the second-order derivatives returned.
  integer, intent(in) :: n
  INTEGER :: I
  REAL, INTENT (IN) , DIMENSION (:)  :: XI, FI
  REAL, INTENT (OUT), DIMENSION (:)  :: P2
  REAL, DIMENSION (n)   :: G, H
  REAL, DIMENSION (n-1) :: D, B, C

  ! Assign the intervals and function differences
  DO I = 1, N
    H(I) = XI(I+1) - XI(I)
    G(I) = FI(I+1) - FI(I)
  END DO

  ! Evaluate the coefficient matrix elements
  DO I = 1, N-1
    D(I) = 2*(H(I+1)+H(I))
    B(I) = 6*(G(I+1)/H(I+1)-G(I)/H(I))
    C(I) = H(I+1)
  END DO

  ! Obtain the second-order derivatives
  CALL TRIDIAGONAL_LINEAR_EQ (N-1, D, C, C, B, G)

  P2(1) = 0
  P2(N+1) = 0
  DO I = 2, N 
    P2(I) = G(I-1)
  END DO

END SUBROUTINE CUBIC_SPLINE


pure SUBROUTINE TRIDIAGONAL_LINEAR_EQ (L, D, E, C, B, Z)
! Function to solve the tridiagonal linear equation set.

  INTEGER, INTENT (IN) :: L
  INTEGER :: I
  REAL, INTENT (IN), DIMENSION (L):: D, E, C, B
  REAL, INTENT (OUT), DIMENSION (L):: Z
  REAL, DIMENSION (L)  :: Y, W
  REAL, DIMENSION (L-1):: V, T

  ! Evaluate the elements in the LU decomposition
  W(1) = D(1)
  V(1) = C(1)
  T(1) = E(1)/W(1)
  DO I = 2, L - 1
    W(I) = D(I)-V(I-1)*T(I-1)
    V(I) = C(I)
    T(I) = E(I)/W(I)
  END DO
  W(L) = D(L)-V(L-1)*T(L-1)

  ! Forward substitution to obtain y
  Y(1) = B(1)/W(1)
  DO I = 2, L
    Y(I) = (B(I)-V(I-1)*Y(I-1))/W(I)
  END DO

  ! Backward substitution to obtain z
  Z(L) = Y(L)
  DO I = L-1, 1, -1
    Z(I) = Y(I) - T(I)*Z(I+1)
  END DO

END SUBROUTINE TRIDIAGONAL_LINEAR_EQ


subroutine usage_stop
  print *,'Usage: cdf2fst -i NetCDF_InputFile.cdf -o FST_OutputFile.fst'
  stop
end subroutine usage_stop


end program cdf2fst