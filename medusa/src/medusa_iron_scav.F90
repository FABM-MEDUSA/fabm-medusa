#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA iron chemsitry and scavenging
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
      logical :: deep_fe_fix
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

   call self%register_implemented_routines((/source_do/))
   call self%get_parameter(self%deep_fe_fix,'deep_fe_fix','stop scavenging for Fe below 0.5 umol / m3 at depths > 1000 m',default=.false.)
   call self%get_parameter(self%xk_FeL,'xk_FeL','(umol m-3)-1','dissociation constant for (Fe+ligand)',default=100._rk)
   call self%get_parameter(self%xLgT,'xLgT','umol m-3','total ligand concentration',default=1._rk)
   call self%get_parameter(self%xk_sc_Fe,'xk_sc_Fe','d-1','scavenging rate of "free" Fe',default=1.e-3_rk,scale_factor=d_per_s)
   call self%get_parameter(self%jiron,'jiron','-','iron scavenging scheme: 1-Dutkiewicz et al. (2005),2-Moore et al. (2004),3-Moore et al. (2008),4-Galbraith et al. (2010)',default=1)

   call self%register_state_dependency(self%id_ZFER,'FER','mmol Fe/m**3', 'iron nutrient')
   call self%register_dependency(self%id_depth, standard_variables%depth)
   call self%register_dependency(self%id_ffastc_loc,'ffastc_loc','mmol C m-2 s-1','local remineralisation of detritus (C)')
   call self%register_dependency(self%id_ffastca_loc,'ffastca_loc','mmol Ca m-2 s-1','local remineralisation of detritus (Ca)')
   call self%register_dependency(self%id_ffastsi_loc,'ffastsi_loc','mmol Si m-2 s-1','local remineralisation of detritus (Si)')
   call self%register_dependency(self%id_fscal_part,'fscal_part','nmol C cm-2 s-1','carbon in suspended particles')
   call self%register_diagnostic_variable(self%id_ffescav,'SCAVENGE','mmol Fe/m**3','scavenged iron')


   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)

   !------------------------------------------------------------------
   ! Iron chemistry and fractionation following the Parekh et al. (2004) scheme adopted 
   ! by the Met.Office, MEDUSA models total iron but considers "free" and ligand-bound forms
   ! for the purposes of scavenging (only "free" iron can be scavenged)
   !------------------------------------------------------------------
   class(type_medusa_iron_scav), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_

   !Local variables
   real(rk) :: ZFER,depth,ffastc,ffastca,ffastsi
   real(rk) :: xFeT,xb_coef_tmp,xb2M4ac,xLgF,xFel,xFeF,xFree,ffescav,xmaxFeF,fdeltaFe
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
   real(rk), parameter :: s_per_d = 86400._rk
   real(rk) :: fbase_scav,fscal_sink,fscal_part,fscal_scav,fscal_csink,fscal_sisink,fscal_casink
   real(rk) :: xCscav1,xCscav2,xk_org,xORGscav,xk_inorg,xINORGscav

    _LOOP_BEGIN_

  !iron chemistry and fractionation

  _GET_(self%id_ZFER,ZFER)
  _GET_(self%id_depth,depth)

  xFeT = ZFER * 1.e3_rk !total iron concentration (mmolFe/m3 -> umolFe/m3)

  ! calculate fractionation (based on Diat-HadOCC; in turn based on Parekh et al., 2004)
  xb_coef_tmp = self%xk_FeL * (self%xLgT - xFeT) - 1.0_rk
  xb2M4ac = max(((xb_coef_tmp * xb_coef_tmp) + (4.0_rk * self%xk_FeL * self%xLgT)), 0._rk)
  ! "free" ligand concentration
  xLgF = 0.5_rk * (xb_coef_tmp + (xb2M4ac**0.5_rk)) / self%xk_FeL
  ! ligand-bound iron concentration
  xFeL = self%xLgT - xLgF
  ! "free" iron concentration (and convert to mmolFe/m3)
  xFeF = (xFeT - xFeL) * 1.e-3_rk
  xFree = xFeF / (ZFER + tiny(ZFER))
  !
  ! scavenging of iron (multiple schemes); Andrew Yool is only really happy with the first one at the moment 
  ! - the others involve assumptions (sometimes guessed at by him) that are potentially questionable
  !
  if (self%jiron == 1) then
  !----------------------------------------------------------------------------------------------------------------
  ! Scheme 1: Dutkiewicz et al. (2005)
  ! This scheme includes a single scavenging term based solely on a fixed rate and the availablility of "free" iron
  !----------------------------------------------------------------------------------------------------------------
     ffescav = self%xk_sc_Fe * xFeF
     xmaxFeF = min((xFeF * 1.e3_rk), 0.3_rk)           ! = umol/m3
     fdeltaFe = (xFeT - (xFeL + xmaxFeF)) * 1.e-3_rk   ! = mmol/m3

     ffescav     = ffescav + fdeltaFe * d_per_s        ! = mmol/m3/d !assuming time scale of fdeltaFe of 1 day

     if (self%deep_fe_fix) then

       if ((depth.ge.1000._rk).and.(xFeT .lt. 0.5_rk)) then
          ffescav = 0._rk
       endif

     endif

  elseif (self%jiron == 2) then
   !---------------------------------------------------------------------------------------------------------------
   ! Scheme 2: Moore et al. (2004)
   ! This scheme includes a single scavenging term that accounts for both suspended and sinking particles in 
   ! the water column; this term scavenges total iron rather than "free" iron
   !---------------------------------------------------------------------------------------------------------------
     _GET_(self%id_ffastc_loc, ffastc)
     _GET_(self%id_fscal_part, fscal_part)

     fbase_scav = 0.12_rk / 365.25_rk * d_per_s
     fscal_sink = ffastc * 1.e2_rk * s_per_d
     fscal_scav = fbase_scav * min(((fscal_sink + fscal_part) / 0.0066_rk), 4._rk)

     if (xFeT .lt. 0.4_rk) then

        fscal_scav = fscal_scav * (xFeT / 0.4_rk)

     elseif (xFeT .gt. 0.6_rk) then

        fscal_scav = fscal_scav + ((xFeT / 0.6_rk) * (6._rk/1.4_rk))

     endif

     ffescav = fscal_scav * ZFER
  elseif (self%jiron == 3) then
   !--------------------------------------------------------------------------------------------------------------
   ! Scheme 3: Moore et al. (2008)
   ! This scheme includes a single scavenging term that accounts for sinking particles in the water column, 
   ! and includes organic C, biogenic opal, calcium carbonate and dust in this (though the latter is 
   ! ignored here); this term scavenges total iron rather than "free" iron
   !--------------------------------------------------------------------------------------------------------------
     _GET_(self%id_ffastc_loc, ffastc)
     _GET_(self%id_ffastsi_loc, ffastsi)
     _GET_(self%id_ffastca_loc, ffastca)

     fbase_scav = 0.00384_rk * d_per_s
     fscal_csink  = ffastc  * 1.e6_rk * 12.011_rk  * 1.e-4_rk ! xmassc = 12.011_rk
     fscal_sisink = ffastsi * 1.e6_rk * 60.084_rk * 1.e-4_rk  ! xmasssi = 60.084_rk
     fscal_casink = ffastca * 1.e6_rk * 100.086_rk * 1.e-4_rk ! xmassca = 100.086_rk

     fscal_sink = (fscal_csink * 6._rk + fscal_sisink + fscal_casink) / (100._rk * 1.e3_rk * d_per_s)

     fscal_scav = fbase_scav * fscal_sink

     if (xFeT .lt. 0.5_rk) then

        fscal_scav = fscal_scav * (xFeT / 0.5_rk)

     elseif (xFeT .gt. 0.6_rk) then

        fscal_scav = fscal_scav + (xFeT- 0.6_rk) * 0.00904_rk
                     
     else
    ! Do nothing

    endif

    ffescav = fscal_scav * ZFER


  elseif (self%jiron == 4) then
  !------------------------------------------------------
  ! Scheme 4: Galbraith et al. (2010)
  ! This scheme includes two scavenging terms, one for organic, particle-based scavenging, and another for 
  ! inorganic scavenging; both terms scavenge "free" iron only
  !------------------------------------------------------
  !
    _GET_(self%id_ffastc_loc, ffastc)
    xCscav1    = ffastc * s_per_d * 12.011_rk / 100._rk
    xCscav2    = (xCscav1 * 1.e-3_rk)**0.58_rk
    xk_org     = 0.5_rk * d_per_s ! ((g C m/3)^-1) / d
    xORGscav   = xk_org * xCscav2 * xFeF
    xk_inorg   = 1000._rk * d_per_s ! ((nmol Fe / kg)^1.5)
    xINORGscav = (xk_inorg * (xFeF * 1026._rk)**1.5_rk) * 1.e-3_rk

    ffescav = xORGscav + xINORGscav

  else

    ffescav = 0._rk

  end if

  _SET_DIAGNOSTIC_(self%id_ffescav,ffescav * s_per_d)
  _SET_ODE_(self%id_ZFER, - ffescav)
   _LOOP_END_

   end subroutine do

  end module medusa_iron_scav
