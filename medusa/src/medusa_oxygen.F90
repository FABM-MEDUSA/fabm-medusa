#include "fabm_driver.h"

!
!************************************************************
!               FABM-MEDUSA oxygen air-sea flux
!************************************************************

module medusa_oxygen

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_oxygen

      ! Variable identifiers
      type (type_state_variable_id)                 :: id_ZOXY
      type (type_horizontal_diagnostic_variable_id) ::  id_fair
      type (type_dependency_id)                     :: id_temp,id_salt
      type (type_horizontal_dependency_id)          :: id_apress,id_wnd

   contains

      procedure :: initialize
      procedure :: do_surface

  end type type_medusa_oxygen

contains

   subroutine initialize(self,configunit)

   class(type_medusa_oxygen),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk 
   call self%register_state_dependency(self%id_ZOXY,'ZOXY','mmol O_2/m**3', 'dissolved oxygen')
   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_apress, standard_variables%surface_air_pressure)
   call self%register_dependency(self%id_wnd,standard_variables%wind_speed)
   call self%register_diagnostic_variable(self%id_fair,'fair','mmol O_2/m^2/d','Air-sea flux of oxygen')

   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

   class(type_medusa_oxygen),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: pt,ps,o2,o2_sato,o2_schmidt,kwo2,o2_sat,pp0,kw660,wnd,o2flux,o2sat
   real(rk) :: a0 = 2.00907_rk
   real(rk) :: a1 = 3.22014_rk
   real(rk) :: a2 = 4.05010_rk
   real(rk) :: a3 = 4.94457_rk
   real(rk) :: a4 = -2.56847E-1_rk
   real(rk) :: a5 = 3.88767_rk
   real(rk) :: b0 = -6.24523E-3_rk
   real(rk) :: b1 = -7.37614E-3_rk
   real(rk) :: b2 = -1.03410E-2_rk
   real(rk) :: b3 = -8.17083E-3_rk
   real(rk) :: c0 = -4.88682E-7_rk
   real(rk) :: tt,tk,ts,ts2,ts3,ts4,ts5
   real(rk) :: ans1, ans2
  ! Wanninkhof (2014) coefficients
   real(rk) :: as0 = 1920.4_rk
   real(rk) :: as1 = -135.6_rk
   real(rk) :: as2 = 5.2121_rk
   real(rk) :: as3 = -0.10939_rk
   real(rk) :: as4 = 0.00093777_rk
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk) :: a(7)
   real(rk) :: b(7)
   real(rk) :: tmp_k
   data a(1) / 0.166_rk /  ! Liss & Merlivat (1986)    [approximated]
   data a(2) / 0.3_rk /    ! Wanninkhof (1992)         [sans enhancement]
   data a(3) / 0.23_rk /   ! Nightingale et al. (2000) [good]
   data a(4) / 0.23_rk /   ! Nightingale et al. (2000) [better]
   data a(5) / 0.222_rk /  ! Nightingale et al. (2000) [best]
   data a(6) / 0.337_rk /  ! OCMIP-2                   [sans variability]
   data a(7) / 0.251_rk /  ! Wanninkhof (2014)         [assumes 6h avg winds]
   data b(1) / 0.133_rk /
   data b(2) / 0.0_rk /
   data b(3) / 0.0_rk /
   data b(4) / 0.1_rk /
   data b(5) / 0.333_rk /
   data b(6) / 0.0_rk /
   data b(7) / 0.0_rk /


   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,pt)
   _GET_(self%id_salt,ps)
   _GET_(self%id_ZOXY,o2)
   _GET_HORIZONTAL_(self%id_apress,pp0)
   _GET_HORIZONTAL_(self%id_wnd,wnd)

      o2 = o2/1000._rk

! Calculate gas transfer velocity (cm/h)

      tmp_k = (a(7) * wnd**2) + (b(7) * wnd)

! Convert tmp_k from cm/h to m/s
      kw660 = tmp_k / (3600._rk * 100._rk)

 !note: air-sea fluxes must be corrected for sea ice

      tt   = 298.15_rk - pt
      tk   = 273.15_rk + pt
      ts   = log(tt / tk)
      ts2  = ts**2_rk
      ts3  = ts**3_rk
      ts4  = ts**4_rk
      ts5  = ts**5_rk

      ans1 = a0 + a1*ts + a2*ts2 + a3*ts3 + a4*ts4 + a5*ts5  &
             + ps*(b0 + b1*ts + b2*ts2 + b3*ts3)             &
             + c0*(ps*ps)
      ans2 = exp(ans1)

!  Convert from ml/l to mol/m3
   o2_sato = (ans2 / 22391.6_rk) * 1000.0_rk

   o2_schmidt = as0 + pt*(as1 + pt*(as2 + pt*(as3 + pt*as4)))
   kwo2 = kw660 * (660._rk / o2_schmidt)**0.5_rk
   o2sat = o2_sato * pp0 / 101325._rk
   !print*,pp0

   o2flux = kwo2 * (o2sat - o2)
   o2flux = o2flux *1000._rk

   _SET_SURFACE_EXCHANGE_(self%id_ZOXY, o2flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fair, o2flux)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module medusa_oxygen
