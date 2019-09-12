#include "fabm_driver.h"

!
!******************************************
!            FABM-MEDUSA optics
!******************************************

module medusa_light

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_light
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZCHN,id_ZCHD
      type (type_dependency_id)            :: id_dz,id_depth
      type (type_diagnostic_variable_id)   :: id_xpar
      type (type_horizontal_diagnostic_variable_id) :: id_qsr_diag, id_MED_XZE
      type (type_horizontal_dependency_id) :: id_qsr
      ! Parameters

      real(rk) :: rpig,xkg0,xkr0,xkgp,xkrp,xlg,xlr

   contains
      procedure :: initialize
      procedure :: get_light
  end type type_medusa_light

contains

   subroutine initialize(self,configunit)

   class(type_medusa_light),intent(inout),target :: self
   integer,                         intent(in)           :: configunit
   call self%get_parameter(self%xkg0,'xkg0','m-1','green water absorption coefficient',default=0.0232_rk)
   call self%get_parameter(self%xkr0,'xkr0','m-1','red water absorption coefficient',default=0.225_rk)
   call self%get_parameter(self%xkgp,'xkgp','m-1','pigment green absorption coefficient',default=0.074_rk)
   call self%get_parameter(self%xkrp,'xkrp','m-1','pigment red absorption coefficient',default=0.037_rk)
   call self%get_parameter(self%xlg,'xlg','-','green chl exposant',default=0.629_rk)
   call self%get_parameter(self%xlr,'xlr','-','red chl exposant',default=0.674_rk)
   call self%get_parameter(self%rpig,'rpig','-','chla / (chla+phea) ratio',default=0.7_rk)

   call self%register_state_dependency(self%id_ZCHN,'CHN','mg chl/m**3', 'chlorophyll in non-diatoms')
   call self%register_state_dependency(self%id_ZCHD,'CHD','mg chl/m**3', 'chlorophyll in diatoms')

   call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
   call self%register_dependency(self%id_qsr, standard_variables%surface_downwelling_shortwave_flux)
   call self%register_diagnostic_variable(self%id_qsr_diag,'MED_QSR','W/m^2','Sea surface radiation',source=source_do_column)
   call self%register_diagnostic_variable(self%id_MED_XZE,'MED_XZE','m','Euphotic depth',source=source_do_column)
   call self%register_diagnostic_variable(self%id_xpar,'MED_XPAR','W/m^2','photosynthetically active radiation', &
        standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux,source=source_do_column)

   end subroutine initialize

   subroutine get_light(self,_ARGUMENTS_VERTICAL_)
   class(type_medusa_light),intent(in) :: self
   
   _DECLARE_ARGUMENTS_VERTICAL_

    real(rk) :: dz,ZCHN,ZCHD,qsr,depth
    integer  :: check
    real(rk) :: totchl,zpar0m,zpar100
    real(rk) :: zpig,zkr,zkg,zparr,zparg,xpar,zparr1,zparg1

    _GET_HORIZONTAL_(self%id_qsr,qsr)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_qsr_diag,qsr)
  
     zpar0m = qsr * 0.43_rk
    
    if (zpar0m .le. 0.0_rk) zpar0m = 0.001_rk

    zpar100 = zpar0m * 0.01_rk

    xpar = zpar0m
    zparr = 0.5_rk * zpar0m
    zparg = 0.5_rk * zpar0m

    depth = 0
    check = 0

   _VERTICAL_LOOP_BEGIN_
    
    _GET_(self%id_dz,dz)
    _GET_(self%id_ZCHN,ZCHN)
    _GET_(self%id_ZCHD,ZCHD)
   
    totchl = ZCHN + ZCHD
    zpig = MAX( TINY(0._rk), totchl/self%rpig) ! total pigment
    zkr  = self%xkr0 + self%xkrp * EXP( self%xlr * LOG( zpig ) ) ! total absorption coefficient in red
    zkg  = self%xkg0 + self%xkgp * EXP( self%xlg * LOG( zpig ) ) ! total absorption coefficient in green
    zparr1 = zparr / zkr / dz * ( 1._rk - EXP( -zkr*dz ) )    ! red compound of par
    zparg1    = zparg / zkg / dz * ( 1._rk - EXP( -zkg*dz ) ) ! green compound of par
    xpar = MAX( zparr1 + zparg1, 1.e-15_rk )
    _SET_DIAGNOSTIC_(self%id_xpar,xpar)
    if ((xpar .lt. zpar100).and.(check==0)) then
        check = 1
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_MED_XZE, depth)
    end if
    zparr = zparr * EXP( -zkr * dz )
    zparg = zparg * EXP( -zkg * dz )

    depth = depth + dz
   _VERTICAL_LOOP_END_

     if (check == 0) _SET_HORIZONTAL_DIAGNOSTIC_(self%id_MED_XZE, depth)

   end subroutine get_light

end module medusa_light
