#include "fabm_driver.h"

!
!********************************************    
!                FABM-MEDUSA                       
!********************************************                                                                            

module fabm_medusa

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZCHN,id_ZCHD,id_ZPHN,id_ZPHD,id_ZDIN,id_ZFER,id_ZSIL
      type (type_dependency_id)            :: id_temp 
  !    type (type_diagnostic_variable_id)   ::

      ! Parameters
      logical :: jliebig
      real(rk) :: xnln,xfln,xnld,xsld,xfld

   contains
      procedure :: initialize
      procedure :: do

  end type

contains

   subroutine initialize(self,configunit)
! SELF%DT???
   class(type_medusa),intent(inout),target :: self
   integer,               intent(in)           :: configunit   
   !Register parameters
   call self%get_parameter(self%xxi, 'xxi', 'molN gC-1','C : N conversion factor', default=0.01257_rk)
   call self%get_parameter(self%jliebig, 'jliebig', 'Liebig''s minimum law for nutrient limitation', default=.false.)
   call self%get_parameter(self%xnln, 'xnln', 'mmolN m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.5_rk)
   call self%get_parameter(self%xfln, 'xfln', 'mmolFe m-3','Fe nutrient uptake half-saturation constant (non-diatoms)', default=0.33_rk)
   call self%get_parameter(self%xnld, 'xnld', 'mmolN m-3','N nutrient uptake half-saturation constant (non-diatoms)', default=0.75_rk)
   call self%get_parameter(self%xsld, 'xsld', 'mmolSi m-3','Si nutrient uptake half-saturation constant (diatoms)', default=3.0_rk)
   call self%get_parameter(self%xfld, 'xfld', 'mmolFe m-3','Fe nutrient uptake half-saturation constant (diatoms)', default=0.67_rk)
   call self%get_parameter(self%xvpn, 'xvpn', 'd-1','Maximum phytoplankton growth rate (non-diatoms)', default=0.53_rk)
   call self%get_parameter(self%xvpd, 'xvpd', 'd-1','Maximum phytoplankton growth rate (diatoms)', default=0.50_rk)
   ! Register state variables
   call self%register_state_variable(self%id_ZCHN,'ZCHN','mg chl/m**3', 'chlorophyll in non-diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZCHD','mg chl/m**3', 'chlorophyll in diatoms', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZPHN','mmolN/m**3', 'non-diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZCHD,'ZPHD','mmolN/m**3', 'diatom phytoplankton', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZDIN,'ZDIN','mmolN/m**3', 'nitrogen nutrient', minimum=0.0_rk) !no river dilution?
   call self%register_state_variable(self%id_ZFER,'ZFER','mmolFe/m**3', 'iron nutrient', minimum=0.0_rk)
   call self%register_state_variable(self%id_ZSIL,'ZSIL','mmolSi/m**3', 'silicic acid', minimum=0.0_rk)

   ! Register diagnostic variables
     
   ! Register environmental dependencies

   call self%register_dependency(self%id_temp, standard_variables%temperature)

   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   class(type_medusa), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

! !LOCAL VARIABLES:
    real(rk) :: fthetan,fthetad !scaled chl/biomass ratio
    real(rk) :: fnln,ffln ! non-diatom Qn/Qf terms
    real(rk) :: fnld,fsld,ffld ! diatom Qn/Qs/Qf terms
    real(rk) :: fun_T,xvpnT,xvpdT ! temperature limitation
    real(rk) :: fpnlim,fpdlim !nutrient limitation of primary production 
    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values

    _GET_(self%id_ZCHN,ZCHN)
    _GET_(self%id_ZCHD,ZCHD)
    _GET_(self%id_ZPHN,ZPHN)
    _GET_(self%id_ZPHD,ZPHD)
    _GET_(self%id_ZDIN,ZDIN)
    _GET_(self%id_ZFER,ZFER)
    _GET_(self%id_ZSIL,ZSIL) 
    _GET_(self%id_temp,temp)
 
   !PHYTOPLANKTON
   !Chlorophyll
   fthetan = (ZCHN * self%xxi) / ZPHN
   fthetad = (ZCHD * self%xxi) / ZPHD
 
  !Temperature limitation
   fun_T = 1.066_rk**temp

   xvpnT = self%xvpn * fun_T
   xvpdT = self%xvpd * fun_T
  
   ! Phytoplankton nutrient limitation
   !! Non-diatoms 
   fnln = ZDIN / (ZDIN + self%xnln) !non-diatom Qn term
   ffln = ZFER / (ZFER + self%xfln) !non-diatom Qf term
   !! Diatoms 
   fnld = ZDIN / (ZDIN + self%xnld) !diatom Qn term
   fsld = ZSIL / (ZSIL + self%xsld) !diatom Qs term
   ffld = ZFER / (ZFER + self%xfld) !diatom Qf term

   ! Primary production (non-diatoms)
   if (self%njliebig=.false.) then
     fpnlim = fnln * ffln
   else
     fpnlim = min(fnln,ffln)
   end if

   !Primary production (diatoms)
   if (self%nliebig=.false.) then
     fpdlim = fnld * fsld
   else
     fpdlim = min(fnld,ffld)
   end if

  !_SET_ODE_(self%id_..,)

   _LOOP_END_

   end subroutine do


  end module fabm_medusa

