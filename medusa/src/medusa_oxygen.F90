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
      type (type_horizontal_diagnostic_variable_id) :: id_fairo2,id_O2SAT
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
   call self%register_dependency(self%id_fr_i,standard_variables%ice_area_fraction)
   ! call self%register_dependency(self%id_apress, standard_variables%surface_air_pressure)
   call self%register_diagnostic_variable(self%id_fairo2,'O2FLUX','mmol O_2/m^2/d','Air-sea O2 flux',source=source_do_surface)
   call self%register_diagnostic_variable(self%id_O2SAT,'O2SAT','mmol O_2/m2/d','Surface ocean O2 saturation',source=source_do_surface)
   call self%register_horizontal_dependency(self%id_kw660, 'KW660', 'm/s', 'gas transfer velocity')

   end subroutine initialize

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

  ! Calculates O2 saturation at 1 atm pressure
  ! Original author : Andrew Yool (14/10/04, revised 08/07/11)
  ! 
  ! This subroutine calculates the oxygen saturation concentration at 1 atmosphere pressure in mol/m3 
  ! given ambient temperature and salinity.  This formulation is (ostensibly) taken from 
  ! Garcia & Gordon (1992, L&O, 1307-1312).  The function works in the range -1.9 <= T <= 40, 0 <= S <= 42.  
  !
  ! Check value : T = 10, S = 35, oxy_sato = 0.282015 mol/m3
  !
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


   ! Convert from ml/l to mol/m3
   o2_sato = (ans2 / 22391.6_rk) * 1000.0_rk

   ! Calculate Schmidt number for ocean uptake of O2
   ! Original author : Andrew Yool (14/10/04, revised 08/07/11)
   ! 
   ! This subroutine calculates the Schmidt number for O2 using sea surface temperature.
   ! The code is based upon that developed as part of the OCMIP-2 project (1998-2000).
   ! The coefficients used are taken from Wanninkhof (2014) 
   ! Winninkhof, R. (2014). Relationship between wind speed and gas
   ! exchange over the ocean revisited. LIMNOLOGY AND OCEANOGRAPHY-METHODS
   ! 12, 351-362, doi:10.4319/lom.2014.12.351
   !
   o2_schmidt = as0 + pt*(as1 + pt*(as2 + pt*(as3 + pt*as4)))
   kwo2 = kw660 * (660._rk / o2_schmidt)**0.5_rk  ! Calculate the transfer velocity for O2 (m/s)

   o2sat = o2_sato * 1._rk ! Calculate the saturation concentration for oxygen (mol/m3)

  ! Calculate time rate of change of O2 due to gas exchange (mol/m3/s), correct for sea-ice
   o2flux = (1._rk - fr_i) * kwo2 * (o2sat - o2)

   o2flux = o2flux *1000._rk      ! Convert oxygen flux and saturation to mmol / m3

   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_O2SAT,o2sat * 1000._rk)   
   _SET_SURFACE_EXCHANGE_(self%id_ZOXY, o2flux)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fairo2, o2flux / d_per_s)

   _HORIZONTAL_LOOP_END_

   end subroutine do_surface

end module medusa_oxygen
