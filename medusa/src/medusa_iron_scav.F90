#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA iron scavenging module
!*********************************************************

module medusa_iron_scav

   use fabm_types

   implicit none

   private

  type,extends(type_base_model),public :: type_medusa_iron_scav
      ! Variable identifiers
      type (type_state_variable_id)        :: id_ZFER
      type (type_dependency_id)            :: id_depth
      type (type_dependency_id)            :: id_ffastc_loc,id_ffastca_loc,id_ffastsi_loc,id_fscal_part
      type (type_diagnostic_variable_id)   :: id_ffescav
  ! Parameters

      real(rk) :: xk_FeL,xLgT,xk_sc_Fe
      integer :: jiron
   contains

      procedure :: initialize
      procedure :: do

  end type type_medusa_iron_scav

contains

   subroutine initialize(self,configunit)

   class(type_medusa_iron_scav),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   call self%get_parameter(self%xk_FeL,'xk_FeL','-','dissociation constant for (Fe+ligand)',default=100._rk)
   call self%get_parameter(self%xLgT,'xLgT','umol m-3','total ligand concentration',default=1._rk)
   call self%get_parameter(self%xk_sc_Fe,'xk_sc_Fe','d-1','scavenging rate of "free" Fe',default=1.E-3_rk,scale_factor=d_per_s)
   call self%get_parameter(self%jiron,'jiron','-','iron scavenging scheme: 1-Dutkiewicz et al. (2005),2-Moore et al. (2004),3-Moore et al. (2008),4-Galbraith et al. (2010)',default=1)

   call self%register_state_dependency(self%id_ZFER,'ZFER','mmol Fe/m**3', 'iron nutrient')
   call self%register_dependency(self%id_depth, standard_variables%depth)
   call self%register_dependency(self%id_ffastc_loc,'ffastc_loc','mmol C m-2 s-1','local remineralisation of detritus (C)')
call self%register_dependency(self%id_ffastca_loc,'ffastca_loc','mmol Ca m-2 s-1','local remineralisation of detritus (Ca)')
call self%register_dependency(self%id_ffastsi_loc,'ffastsi_loc','mmol Si m-2 s-1','local remineralisation of detritus (Si)')
   call self%register_dependency(self%id_fscal_part,'fscal_part','nmol C cm-2 s-1','carbon in suspended particles')
   call self%register_diagnostic_variable(self%id_ffescav,'scav_flux','mmol Fe/m**3','scavenged iron')


   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   class(type_medusa_iron_scav), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

   !Local variables
   real(rk) :: ZFER,depth,ffastc,ffastca,ffastsi
   real(rk) :: xFeT,xb_coef_tmp,xb2M4ac,xLgF,xFel,xFeF,xFree,ffescav,xmaxFeF,fdeltaFe
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk) :: fbase_scav,fscal_sink,fscal_part,fscal_scav,fscal_csink,fscal_sisink,fscal_casink
   real(rk) :: xCscav1,xCscav2,xk_org,xORGscav,xk_inorg,xINORGscav

    _LOOP_BEGIN_

  !iron chemistry and fractionation

  _GET_(self%id_ZFER,ZFER)
  _GET_(self%id_depth,depth)

  xFeT = ZFER * 1.e3_rk !total iron concentration (mmolFe/m3 -> umolFe/m3)
  xb_coef_tmp = self%xk_FeL * (self%xLgT - xFeT) - 1.0_rk
  xb2M4ac = max(((xb_coef_tmp * xb_coef_tmp) + (4.0_rk * self%xk_FeL * self%xLgT)), 0._rk)
  xLgF = 0.5_rk * (xb_coef_tmp + (xb2M4ac**0.5_rk)) / self%xk_FeL ! "free" ligand concentration
  xFeL = self%xLgT - xLgF ! ligand-bound iron concentration
  xFeF = (xFeT - xFeL) * 1.e-3_rk! "free" iron concentration (and convert to mmolFe/m3)
  xFree = xFeF / (ZFER + tiny(ZFER))

 ! ! Scavenging of iron
  if (self%jiron == 1) then
     ffescav = self%xk_sc_Fe * xFeF
     xmaxFeF = min((xFeF * 1.e3_rk), 0.3_rk)        ! = umol/m3
     fdeltaFe = (xFeT - (xFeL + xmaxFeF)) * 1.e-3_rk   ! = mmol/m3

     ffescav     = ffescav + fdeltaFe / d_per_s        ! = mmol/m3/d !assuming time scale of fdeltaFe of 1 day

   !  if ((depth .gt. 1000._rk) .and. (xFeT .lt. 0.5_rk)) then

   !     ffescav = 0._rk

   !  endif

  elseif (self%jiron == 2) then

     _GET_(self%id_ffastc_loc, ffastc)
     _GET_(self%id_fscal_part, fscal_part)

     fbase_scav = 0.12_rk / 365.5_rk * d_per_s
     fscal_sink = ffastc * 1.e2_rk
     fscal_scav = fbase_scav * min(((fscal_sink + fscal_part) / 0.0066_rk), 4._rk)
     if (xFeT .lt. 0.4_rk) then

        fscal_scav = fscal_scav * (xFeT / 0.4_rk)

     elseif (xFeT .gt. 0.6_rk) then

        fscal_scav = fscal_scav + ((xFeT / 0.6_rk) * (6._rk/1.4_rk))

     endif

     ffescav = fscal_scav * ZFER
  elseif (self%jiron == 3) then

     _GET_(self%id_ffastc_loc, ffastc)
     _GET_(self%id_ffastsi_loc, ffastsi)
     _GET_(self%id_ffastca_loc, ffastca)

     fbase_scav = 0.00384_rk
     fscal_csink  = ffastc  * 1.e6_rk * 12.011_rk  * 1.e-4_rk * d_per_s ! Using 12.011 rather than refer to parameter xmassc
     fscal_sisink = ffastsi * 1.e6_rk * 60.084_rk * 1.e-4_rk * d_per_s !xmasssi = 60.084_rk
     fscal_casink = ffastca * 1.e6_rk * 100.086_rk * 1.e-4_rk * d_per_s !xmassca = 100.086_rk

     fscal_sink = (fscal_csink * 6._rk + fscal_sisink + fscal_casink) / (100._rk * 1.e3_rk * d_per_s)

     if (xFeT .lt. 0.5_rk) then

        fscal_scav = fscal_scav * xFeT / 0.5_rk

     elseif (xFeT .gt. 0.6_rk) then

        fscal_scav = fscal_scav + (xFeT- 0.6_rk) * 0.00904_rk
                     
     else
    ! Do nothing... what do you mean do nothing??

    endif

    ffescav = fscal_scav * ZFER


  elseif (self%jiron == 4) then

    _GET_(self%id_ffastc_loc, ffastc)
    xCscav1    = (ffastc * 12.011_rk) / (100._rk * d_per_s)
    xCscav2    = (xCscav1 * 1.e-3_rk)**0.58_rk
    xk_org     = 0.5_rk * d_per_s ! ((g C m/3)^-1) / d
    xORGscav   = xk_org * xCscav2 * xFeF
    xk_inorg   = 1000._rk * d_per_s ! ((nmol Fe / kg)^1.5)
    xINORGscav = (xk_inorg * (xFeF * 1026._rk)**1.5_rk) * 1.e-3_rk

    ffescav = xORGscav + xINORGscav

  else

    ffescav = 0._rk

  end if

  _SET_DIAGNOSTIC_(self%id_ffescav,ffescav*86400._rk)
  _SET_ODE_(self%id_ZFER, - ffescav)
   _LOOP_END_

   end subroutine do

  end module medusa_iron_scav
