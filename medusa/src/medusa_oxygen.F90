#include "fabm_driver.h"

!
!************************************************************
!               FABM-MEDUSA oxygen air-sea flux
!************************************************************

module medusa_oxygen

   use fabm_types
   use fabm_standard_variables

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_oxygen

      ! Variable identifiers
      type (type_state_variable_id)                 :: id_ZOXY
      type (type_horizontal_diagnostic_variable_id) :: id_fairo2
      type (type_dependency_id)                     :: id_temp,id_salt
      type (type_horizontal_dependency_id)          :: id_kw660, id_fr_i !id_apress
          
   contains

      procedure :: initialize
      procedure :: do_surface

  end type type_medusa_oxygen

contains

   subroutine initialize(self,configunit)

   class(type_medusa_oxygen),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk 
   call self%register_state_dependency(self%id_ZOXY,'OXY','mmol O_2/m**3', 'dissolved oxygen')
   ! Register environmental dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_salt, standard_variables%practical_salinity)
   call self%register_dependency(self%id_fr_i, type_horizontal_standard_variable(name='ice_fraction'))
   ! call self%register_dependency(self%id_apress, standard_variables%surface_air_pressure)
   call self%register_diagnostic_variable(self%id_fairo2,'fairo2','mmol O_2/m^2/d','Air-sea flux of oxygen')

   call self%register_horizontal_dependency(self%id_kw660, 'kw660', 'm/s', 'gas transfer velocity')

   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

   class(type_medusa_oxygen),intent(in) :: self

   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: pt,ps,o2,o2_sato,o2_schmidt,kwo2,o2_sat,pp0,kw660,o2flux,o2sat,fr_i
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

   _HORIZONTAL_LOOP_BEGIN_

   _GET_(self%id_temp,pt)
   _GET_(self%id_salt,ps)
   _GET_(self%id_ZOXY,o2)
   !_GET_HORIZONTAL_(self%id_apress,pp0)
   _GET_HORIZONTAL_(self%id_kw660,kw660)
   _GET_HORIZONTAL_(self%id_fr_i,fr_i)
      o2 = o2/1000._rk

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

   o2sat = o2_sato * 1._rk !Use this value for now !* pp0 / 101325._rk

   o2flux = (1._rk - fr_i) * kwo2 * (o2sat - o2)
   o2flux = o2flux *1000._rk
   
   _SET_SURFACE_EXCHANGE_(self%id_ZOXY, o2flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fairo2, o2flux*86400._rk)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module medusa_oxygen
