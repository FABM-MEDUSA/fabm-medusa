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
      type (type_state_variable_id)               :: id_ZDIN,id_ZSIL,id_ZFER,id_ZDIC,id_ZDET,id_ZDTC,id_ZOXY,id_ZALK
     ! type (type_dependency_id)            ::
     ! type (type_diagnostic_variable_id)   ::

      ! Parameters
      real(rk) :: xsedn,xsedc,xsedfe,xsedsi,xsedca
      real(rk) :: xrfn
      real(rk) :: xthetanit,xthetarem,xo2min

   contains

      procedure :: initialize
      procedure :: do_bottom

  end type type_medusa_benthic

contains

   subroutine initialize(self,configunit)

   class(type_medusa_benthic),intent(inout),target :: self
   integer,               intent(in)           :: configunit

   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
   call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
   call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)
   call self%get_parameter(self%xsedc, 'xsedc', 'd-1','benthic C remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedn, 'xsedn', 'd-1','benthic N remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedfe, 'xsedfe', 'd-1','benthic Fe remineralisation rate', default=0.05_rk,scale_factor=d_per_s) !NB: check default value
   call self%get_parameter(self%xsedsi, 'xsedsi', 'd-1','benthic Si remineralisation rate', default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xsedca, 'xsedca', 'd-1','benthic CaCO3 remineralisation rate', default=0.01_rk,scale_factor=d_per_s)
   call self%get_parameter(self%xrfn,'xrfn','umol Fe mol N-1 m','phytoplankton Fe : N uptake ratio',default=30.0e-6_rk)

   call self%register_state_variable(self%id_ZSEDC,'BEN_C','mmol C/m**2', 'sediment organic carbon pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDN,'BEN_N','mmol N/m**2', 'sediment organic nitrogen pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDFE,'BEN_FE','mmol Fe/m**2', 'sediment organic iron pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDSI,'BEN_SI','mmol Si/m**2', 'sediment organic silica pool', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSEDCA,'BEN_CA','mmol C/m**2', 'sediment organic calcite pool', minimum=0.0_rk)

   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_ZSEDC)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_ZSEDN)
   call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_ZSEDSI)

   call self%register_state_dependency(self%id_ZOXY,'OXY','mmol O_2 m-3', 'dissolved oxygen')
   call self%register_state_dependency(self%id_ZDIN,'DIN','mmol N m-3', 'nitrogen nutrient')
   call self%register_state_dependency(self%id_ZSIL,'SIL','mmol Si m-3', 'silicic acid')
   call self%register_state_dependency(self%id_ZFER,'FER','mmol Fe m-3', 'iron nutrient')
   call self%register_state_dependency(self%id_ZDIC,'DiC','mmol C m-3', 'dissolved inorganic carbon')
   call self%register_state_dependency(self%id_ZDET,'DET','mmol N m-3', 'detritus nitrogen')
   call self%register_state_dependency(self%id_ZDTC,'DTC','mmol C m-3', 'detritus carbon')
   call self%register_state_dependency(self%id_ZALK,'ALK','meq m-3', 'total alkalinity')

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_)

   class(type_medusa_benthic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_

! !LOCAL VARIABLES:

    real(rk) :: ZSEDC,ZSEDN,ZSEDFE,ZSEDSI,ZSEDCA,ZDET,ZDTC,ZSIL,ZALK,ZOXY
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _HORIZONTAL_LOOP_BEGIN_

    _GET_HORIZONTAL_(self%id_ZSEDC,ZSEDC)
    _GET_HORIZONTAL_(self%id_ZSEDN,ZSEDN)
    _GET_HORIZONTAL_(self%id_ZSEDFE,ZSEDFE)
    _GET_HORIZONTAL_(self%id_ZSEDSI,ZSEDSI)
    _GET_HORIZONTAL_(self%id_ZSEDCA,ZSEDCA)
    _GET_(self%id_ZDET,ZDET)
    _GET_(self%id_ZDTC,ZDTC)
    _GET_(self%id_ZOXY,ZOXY)

    _SET_BOTTOM_ODE_(self%id_ZSEDC,  -self%xsedc * ZSEDC)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, + self%xsedc * ZSEDC)

    _SET_BOTTOM_ODE_(self%id_ZSEDN,  -self%xsedn * ZSEDN)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIN, + self%xsedn * ZSEDN)
     
     if (ZOXY .ge. self%xo2min) _SET_BOTTOM_EXCHANGE_(self%id_ZOXY, - self%xthetanit * self%xsedn * ZSEDN - self%xthetarem * self%xsedc * ZSEDC)

    _SET_BOTTOM_ODE_(self%id_ZSEDFE, -self%xsedfe * ZSEDFE)
    _SET_BOTTOM_EXCHANGE_(self%id_ZFER, + self%xsedfe * ZSEDFE)

    _SET_BOTTOM_ODE_(self%id_ZSEDSI, -self%xsedsi * ZSEDSI)
    _SET_BOTTOM_EXCHANGE_(self%id_ZSIL, + self%xsedsi * ZSEDSI)

    _SET_BOTTOM_ODE_(self%id_ZSEDCA, -self%xsedca * ZSEDCA)
    _SET_BOTTOM_EXCHANGE_(self%id_ZALK, + 2._rk * self%xsedca * ZSEDCA)
    _SET_BOTTOM_EXCHANGE_(self%id_ZDIC, self%xsedca * ZSEDCA)

    _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

  end module medusa_benthic
