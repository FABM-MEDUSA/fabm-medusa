#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA benthic module
!*********************************************************

module medusa_benthic

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_benthic
      ! Variable identifiers
      type (type_bottom_state_variable_id)        :: id_ZSEDC,id_ZSEDN,id_ZSEDFE,id_ZSEDSI,id_ZSEDCA
     ! type (type_dependency_id)            ::
     ! type (type_diagnostic_variable_id)   ::
     ! type (type_state_variable_id)   ::

      ! Parameters
      real(rk) :: xsedn,xsedc,xsedfe,xsedsi,xsedca
      real(rk) :: wdep,xrfn


   contains

      procedure :: initialize
      procedure :: do_bottom

  end type type_medusa_benthic

contains

   subroutine initialize(self,configunit)

   class(type_medusa_benthic),intent(inout),target :: self
   integer,               intent(in)           :: configunit

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   call self%get_parameter(self%xsedc, 'xsedc', 'd-1','benthic C remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedn, 'xsedn', 'd-1','benthic N remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedfe, 'xsedfe', 'd-1','benthic Fe remineralisation rate', default=0.01_rk,scale_factor=d_per_s) !NB: check default value
   call self%get_parameter(self%xsedsi, 'xsedsi', 'd-1','benthic Si remineralisation rate', default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedca, 'xsedca', 'd-1','benthic CaCO3 remineralisation rate', default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%wdep,'wdep','m d-1','detritus deposition rate', default=2.5_rk, scale_factor=d_per_s)
   call self%get_parameter(self%xrfn,'xrfn','umol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=0.03_rk)

   call self%register_state_variable(self%id_ZSEDC,'ZSEDC','mmol C/m**2', 'sediment organic carbon pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDN,'ZSEDN','mmol N/m**2', 'sediment organic nitrogen pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDFE,'ZSEDFE','mmol Fe/m**2', 'sediment organic iron pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDC,'ZSEDSI','mmol Si/m**2', 'sediment organic silica pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDC,'ZSEDCA','mmol C/m**2', 'sediment organic calcite pool', minimum=0.0_rk)

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_)

   class(type_medusa_benthic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_

! !LOCAL VARIABLES:

    real(rk) :: ZSEDC,ZSEDN,ZSEDFE,ZSEDSI,ZSEDCA,ZDET,ZDTC,ZSIL,ZALK
    real(rk) :: fluxc,fluxn,fluxfe
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _HORIZONTAL_LOOP_BEGIN_

    _GET_HORIZONTAL_(self%id_ZSEDC,ZSEDC)
    _GET_HORIZONTAL_(self%id_ZSEDN,ZSEDN)
    _GET_HORIZONTAL_(self%id_ZSEDFE,ZSEDFE)
    _GET_HORIZONTAL_(self%id_ZSEDSI,ZSEDSI)
    _GET_HORIZONTAL_(self%id_ZSEDCA,ZSEDCA)
    _GET_(self%id_ZDET,ZDET)
    _GET_(self%id_ZDTC,ZDTC)
    _GET_(self%id_ZSIL,ZSIL)
    _GET_(self%id_ZALK,ZALK)

     fluxc = self%wdep * ZDTC
     fluxn = self%wdep * ZDET
     fluxfe = self%wdep * ZDET * self%xrfn

    _SET_BOTTOM_EXCHANGE_(self%id_ZDET,-fluxn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDTC,-fluxc)

    _SET_BOTTOM_ODE_(self%id_ZSEDC,  -self%xsedc * ZSEDC + fluxc)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, + self%xsedc * ZSEDC)

    _SET_BOTTOM_ODE_(self%id_ZSEDN,  -self%xsedn * ZSEDN + fluxn)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIN, + self%xsedn * ZSEDN)

    _SET_BOTTOM_ODE_(self%id_ZSEDFE, -self%xsedfe * ZSEDFE + fluxfe)

    _SET_BOTTOM_ODE_(self%id_ZSEDSI, -self%xsedsi * ZSEDSI)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + self%xsedsi * ZSEDSI)

    _SET_BOTTOM_ODE_(self%id_ZSEDCA, -self%xsedca * ZSEDCA)
    _SET_BOTTOM_EXCHANGE_(self%id_ZALK, + self%xsedca * ZSEDCA)

    _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

  end module medusa_benthic
