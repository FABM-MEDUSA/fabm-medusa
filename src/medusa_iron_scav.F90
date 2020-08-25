#include "fabm_driver.h"

!
!*********************************************************
!            FABM-MEDUSA iron chemsitry and scavenging
!*********************************************************

module medusa_iron_scav

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_medusa_iron_scav
      ! Variable identifiers
      type (type_state_variable_id)      :: id_ZFER
      type (type_dependency_id)          :: id_depth
      type (type_dependency_id)          :: id_ffastc_loc,id_ffastca_loc,id_ffastsi_loc,id_fscal_part
      type (type_diagnostic_variable_id) :: id_ffescav

      ! Parameters
      logical :: deep_fe_fix
      real(rk) :: xk_FeL,xLgT,xk_sc_Fe
      integer :: jiron
   contains
      procedure :: initialize
      procedure :: do
   end type type_medusa_iron_scav

contains

   subroutine initialize(self, configunit)
      class(type_medusa_iron_scav), intent(inout), target :: self
      integer,                      intent(in)            :: configunit
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

   subroutine do(self, _ARGUMENTS_DO_)
      class(type_medusa_iron_scav), intent(in) :: self
     _DECLARE_ARGUMENTS_DO_

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk
      real(rk), parameter :: s_per_d = 86400._rk
      real(rk) :: ZFER,depth,ffastc,ffastca,ffastsi
      ! parametes for Parekh et al (2005) iron scheme
      real(rk) :: xFeT,xLgF,xFel,xFeF,xFree,ffescav
      ! iron-ligand parameters
      real(rk) :: xb_coef_tmp,xb2M4ac
      ! max Fe' parameters
      real(rk) :: xmaxFeF,fdeltaFe
      ! local parameters for Moore et al (2004) alternative scavenging scheme
      real(rk) :: fbase_scav,fscal_sink,fscal_part,fscal_scav
      ! local parameters for Moore et al (2008) alternative scavenging scheme
      real(rk) :: fscal_csink,fscal_sisink,fscal_casink
      ! local parameters for Galbraith et al (2010) alternative scavenging scheme
      real(rk) :: xCscav1,xCscav2,xk_org,xORGscav,xk_inorg,xINORGscav

      _LOOP_BEGIN_

         !------------------------------------------------------------------
         ! Iron chemistry and fractionation
         ! following the Parekh et al. (2004) scheme adopted by the Met.
         ! Office, MEDUSA models total iron but considers "free" and
         ! ligand-bound forms for the purposes of scavenging (only "free"
         ! iron can be scavenged)
         !------------------------------------------------------------------

         _GET_(self%id_ZFER,ZFER)
         _GET_(self%id_depth,depth)

         !total iron concentration (mmol Fe / m3 -> umol Fe / m3)
         xFeT = ZFER * 1.e3_rk

         ! calculate fractionation (based on Diat-HadOCC; in turn based on Parekh et al., 2004)
         xb_coef_tmp = self%xk_FeL * (self%xLgT - xFeT) - 1.0_rk
         xb2M4ac = max(((xb_coef_tmp * xb_coef_tmp) + (4.0_rk * self%xk_FeL * self%xLgT)), 0._rk)

         ! "free" ligand concentration
         xLgF = 0.5_rk * (xb_coef_tmp + (xb2M4ac**0.5_rk)) / self%xk_FeL

         ! ligand-bound iron concentration
         xFeL = self%xLgT - xLgF

         ! "free" iron concentration (and convert to mmol Fe / m3)
         xFeF = (xFeT - xFeL) * 1.e-3_rk
         xFree = xFeF / (ZFER + tiny(ZFER))
         !
         ! scavenging of iron (multiple schemes); the original MEDUSA authors
         ! are only really happy with the first one at the moment - the others involve
         ! assumptions (sometimes guessed at) that are potentially questionable
         !
         if (self%jiron == 1) then
            !----------------------------------------------------------------------------
            ! Scheme 1: Dutkiewicz et al. (2005)
            ! This scheme includes a single scavenging term based
            ! solely on a fixed rate and the availablility of "free" iron
            !----------------------------------------------------------------------------
            ffescav = self%xk_sc_Fe * xFeF                    ! = mmol/m3/s

            ! Mick Follows code contains a further (optional) implicit
            ! "scavenging" of iron that sets an upper bound on
            ! "free" iron concentration, and essentially caps the
            ! concentration of total iron as xFeL + "free" iron;
            ! since the former is constrained by a fixed total
            ! ligand concentration (= 1.0 umol/m3), and the latter
            ! isn't allowed above this upper bound, total iron is
            ! constrained to a maximum of ...
            !
            !    xFeL(ji,jj) + min(xFeF(ji,jj), 0.3 umol/m3) = 1.0 + 0.3
            !                                  = 1.3 umol / m3
            !
            ! In Mick's code, the actual value of total iron is
            ! reset to this sum (i.e. TFe = FeL + Fe'; but
            ! Fe' <= 0.3 umol/m3)
            ! However, here the amount scavenged is augmented
            ! by an additional amount that serves to drag total
            ! iron back towards that expected from this limitation
            ! on iron concentration ...
            xmaxFeF = min((xFeF * 1.e3_rk), 0.3_rk)           ! = umol/m3

            ! Here, the difference between current total Fe and
            ! (FeL + Fe') is calculated and added to the scavenging
            ! flux already calculated above ...
            fdeltaFe = (xFeT - (xFeL + xmaxFeF)) * 1.e-3_rk   ! = mmol/m3

            ! This assumes that the "excess" iron is dissipated
            ! with a time-scale of 1 day;
            ffescav     = ffescav + fdeltaFe * d_per_s ! = mmol/m3/s

            if (self%deep_fe_fix) then
               ! stop scavenging for iron concentrations below
               ! 0.5 umol / m3 at depths greater than 1000 m; this
               ! aims to end MEDUSA's continual loss of iron at depth
               ! without impacting things at the surface too much; the
               ! justification for this is that it appears to be what
               ! Mick Follows et al. do in their work; it looks like 
               ! this seemingly arbitrary approach effectively
               ! "parameterises" the particle-based scavenging rates
               ! that other models use (i.e. at depth there are no 
               !sinking particles, so scavenging stops)
               if ((depth.ge.1000._rk).and.(xFeT .lt. 0.5_rk)) then
                  ffescav = 0._rk
               endif
           endif

         elseif (self%jiron == 2) then
            !-------------------------------------------------------------------------
            ! Scheme 2: Moore et al. (2004)
            ! This scheme includes a single scavenging term that
            ! accounts for both suspended and sinking particles in
            ! the water column; this term scavenges total iron rather
            ! than "free" iron
            !-------------------------------------------------------------------------
            _GET_(self%id_ffastc_loc, ffastc)
            _GET_(self%id_fscal_part, fscal_part)

            ! this has a base scavenging rate (12% / y) which is
            ! modified by local particle concentration and sinking
            ! flux and which is accelerated when Fe concentration gets
            ! 0.6 nM (= 0.6 umol/m3 = 0.0006 mmol/m3), and decreased
            ! as concentrations below 0.4 nM (= 0.4 umol/m3 =
            ! 0.0004 mmol/m3)

            ! base scavenging rate (0.12 / y)
            fbase_scav = 0.12_rk / 365.25_rk * d_per_s

            ! calculate sinking particle part of scaling factor
            ! this takes local fast sinking carbon (mmol C / m2 / d)
            ! and gets it into nmol C / cm3 / d
            fscal_sink = ffastc * 1.e2_rk * s_per_d

            ! retrieve particle part of scaling factor fscal_part from
            ! medusa_pelagic. This totals up the carbon in suspended particles
            ! (Pn, Pd, Zmi, Zme, D), which comes out in mmol C / m3 (= nmol C / cm3), 
            ! and then multiplies it by a magic factor, 0.002, to get it
            ! into nmol C / cm2 / s
   
            ! calculate scaling factor for base scavenging rate
            ! this uses the sinking flux and standing particle concentration,
            ! divides through by some sort of reference value (= 0.0066 nmol C / cm2 / s)
            ! and then uses this, or not if its too high, to rescale the
            ! base scavenging rate
            fscal_scav = fbase_scav * min(((fscal_sink + fscal_part) / 0.0066_rk), 4._rk)

            ! the resulting scavenging rate is then scaled further
            ! according to the local iron concentration (i.e.
            ! diminished in low iron regions; enhanced in high iron
            ! regions; less alone in intermediate iron regions)
            if (xFeT .lt. 0.4_rk) then !low iron region

               fscal_scav = fscal_scav * (xFeT / 0.4_rk)

            elseif (xFeT .gt. 0.6_rk) then !high iron region

               fscal_scav = fscal_scav + ((xFeT / 0.6_rk) * (6._rk/1.4_rk))

            endif

            ! apply the calculated scavenging rate
            ffescav = fscal_scav * ZFER

         elseif (self%jiron == 3) then

            !---------------------------------------------------------------------
            ! Scheme 3: Moore et al. (2008)
            ! This scheme includes a single scavenging term that
            ! accounts for sinking particles in the water column, 
            ! and includes organic C, biogenic opal, calcium
            ! carbonate and dust in this (though the latter is 
            ! ignored here); this term scavenges total iron
            ! rather than "free" iron
            !---------------------------------------------------------------------
            _GET_(self%id_ffastc_loc, ffastc)
            _GET_(self%id_ffastsi_loc, ffastsi)
            _GET_(self%id_ffastca_loc, ffastca)

            ! this has a base scavenging rate which is modified by
            ! local particle sinking flux and which is accelerated
            ! when Fe concentration is > 0.6 nM (= 0.6 umol/m3 =
            ! 0.0006 mmol/m3), and decreased as concentrations <
            ! 0.5 nM (= 0.5 umol/m3 = 0.0005 mmol/m3)
            fbase_scav = 0.00384_rk * d_per_s

            ! calculate sinking particle part of scaling factor;
            ! this converts mmol / m2 / s fluxes of organic carbon,
            ! silicon and calcium carbonate into ng / cm2 / s
            ! fluxes; it is assumed here that the mass conversions
            ! simply consider the mass of the main element
            ! (C, Si and Ca) and *not* the mass of the molecules
            ! that they are part of; Moore et al. (2008) is unclear
            ! on the conversion that should be used
            !
            ! milli -> nano; mol -> gram; /m2 -> /cm2;
            ! ng C  / cm2 / s
            fscal_csink  = ffastc  * 1.e6_rk * 12.011_rk  * 1.e-4_rk ! xmassc = 12.011_rk
            ! ng Si / cm2 / s
            fscal_sisink = ffastsi * 1.e6_rk * 60.084_rk * 1.e-4_rk  ! xmasssi = 60.084_rk
            ! ng Ca / cm2 / s
            fscal_casink = ffastca * 1.e6_rk * 100.086_rk * 1.e-4_rk ! xmassca = 100.086_rk

            !sum up these sinking fluxes and convert to ng / cm
            ! by dividing through by a sinking rate of
            ! 100 m / d = 1.157 cm / s
            ! ng / cm
            fscal_sink = (fscal_csink * 6._rk + fscal_sisink + fscal_casink) / (100._rk * 1.e3_rk * d_per_s)
 
            ! now calculate the scavenging rate based upon the base
            ! rate and this particle flux scaling; according to the
            ! published units, the result actually has *no* units,
            ! but as it must be expressed per unit time for it to
            ! make any sense
            fscal_scav = fbase_scav * fscal_sink

            ! the resulting scavenging rate is then scaled further
            ! according to the local iron concentration (i.e.
            ! diminished in low iron regions; enhanced in high iron
            ! regions; less alone in intermediate iron regions)
            if (xFeT .lt. 0.5_rk) then 
            ! low iron region (0.5 instead of the 0.4 in Moore et al., 2004)

               fscal_scav = fscal_scav * (xFeT / 0.5_rk)

            elseif (xFeT .gt. 0.6_rk) then
            ! high iron region (functional form different in Moore et al., 2004)

               fscal_scav = fscal_scav + (xFeT- 0.6_rk) * 0.00904_rk
                     
            else
            ! intermediate iron region; do nothing

            endif

            ! apply the calculated scavenging rate
            ffescav = fscal_scav * ZFER

         elseif (self%jiron == 4) then
            !------------------------------------------------------
            ! Scheme 4: Galbraith et al. (2010)
            ! This scheme includes two scavenging terms, one for
            ! organic, particle-based scavenging, and another for 
            ! inorganic scavenging; both terms scavenge "free" iron only
            !------------------------------------------------------
            !
            _GET_(self%id_ffastc_loc, ffastc)

            ! Galbraith et al. (2010) present a more straightforward
            ! outline of the scheme in Parekh et al. (2005) ...
  
            ! sinking particulate carbon available for scavenging
            ! this assumes a sinking rate of 100 m / d (Moore & Braucher, 2008)
            xCscav1    = ffastc * s_per_d * 12.011_rk / 100._rk ! = mg C / m3 / d

            ! scale by Honeyman et al. (1981) exponent coefficient
            ! multiply by 1.e-3 to express C flux in g C rather than mg C
            xCscav2    = (xCscav1 * 1.e-3_rk)**0.58_rk

            ! multiply by Galbraith et al. (2010) scavenging rate
            xk_org     = 0.5_rk * d_per_s ! ((g C m/3)^-1) / s
            xORGscav   = xk_org * xCscav2 * xFeF

            ! Galbraith et al. (2010) also include an inorganic bit ... 
            ! this occurs at a fixed rate, again based on the
            ! availability of "free" iron
            !
            ! k_inorg = 1000 d**-1 nmol Fe**-0.5 kg**-0.5
            !
            ! to implement this here, scale xFeF by 1026 to put in
            ! units of umol/m3 which approximately equal nmol/kg
            xk_inorg   = 1000._rk * d_per_s ! ((nmol Fe / kg)^1.5)
            xINORGscav = (xk_inorg * (xFeF * 1026._rk)**1.5_rk) * 1.e-3_rk

            ! sum these two terms together
            ffescav = xORGscav + xINORGscav

         else

            ! No scheme
            ffescav = 0._rk

         end if

         _SET_DIAGNOSTIC_(self%id_ffescav,ffescav * s_per_d)
         _ADD_SOURCE_(self%id_ZFER, - ffescav)

      _LOOP_END_

   end subroutine do

  end module medusa_iron_scav
