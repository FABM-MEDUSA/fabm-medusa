#include "fabm_driver.h"

!
!*************************************************************
!            FABM-MEDUSA benthic remineralisation module
!*************************************************************

module medusa_benthic

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_medusa_benthic
      ! Variable identifiers
      type (type_bottom_state_variable_id)      :: id_ZSEDC,id_ZSEDN,id_ZSEDFE,id_ZSEDSI,id_ZSEDCA
      type (type_state_variable_id)             :: id_ZDIN,id_ZSIL,id_ZFER,id_ZDIC,id_ZDET,id_ZDTC,id_ZOXY,id_ZALK
      type (type_bottom_diagnostic_variable_id) :: id_f_benout_c,id_f_benout_n,id_f_benout_fe
      type (type_bottom_diagnostic_variable_id) :: id_f_benout_ca,id_f_benout_si,id_f_lyso_ca
      type (type_horizontal_dependency_id)      :: id_CAL_CCD
      type (type_dependency_id)                 :: id_depth

      ! Parameters
      real(rk) :: xsedn,xsedc,xsedfe,xsedsi,xsedca
      real(rk) :: xrfn
      real(rk) :: xthetanit,xthetarem,xo2min
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type type_medusa_benthic

contains

   subroutine initialize(self, configunit)
      class(type_medusa_benthic), intent(inout), target :: self
      integer,                    intent(in)            :: configunit

      real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

      call self%register_implemented_routines((/source_do_bottom/))
      call self%get_parameter(self%xthetanit,'xthetanit','mol O_2 mol N-1','O2 consumption by N remineralisation',default=2.0_rk)
      call self%get_parameter(self%xthetarem,'xthetarem','mol O_2 mol C-1','O2 consumption by C remineralisation',default=1.1226_rk)
      call self%get_parameter(self%xo2min,'xo2min','mmol O_2 m-3','minimum O2 concentration',default=4.0_rk)
      call self%get_parameter(self%xsedc, 'xsedc', 'd-1','benthic C remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
      call self%get_parameter(self%xsedn, 'xsedn', 'd-1','benthic N remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
      call self%get_parameter(self%xsedfe, 'xsedfe', 'd-1','benthic Fe remineralisation rate', default=0.05_rk,scale_factor=d_per_s)
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

      call self%register_dependency(self%id_CAL_CCD,'CAL_CCD','m','calcite CCD depth')
      call self%register_dependency(self%id_depth, standard_variables%depth)

      call self%register_diagnostic_variable(self%id_f_benout_c,'OBEN_C','mmolC/m2/d','Benthic output carbon')
      call self%register_diagnostic_variable(self%id_f_benout_n,'OBEN_N','mmolN/m2/d','Benthic output nitrogen')
      call self%register_diagnostic_variable(self%id_f_benout_fe,'OBEN_FE','mmolFe/m2/d','Benthic output iron')
      call self%register_diagnostic_variable(self%id_f_benout_si,'OBEN_SI','mmolSi/m2/d','Benthic output silicate')
      call self%register_diagnostic_variable(self%id_f_benout_ca,'OBEN_CA','mmolCa/m2/d','Benthic output CaCO3')
      call self%register_diagnostic_variable(self%id_f_lyso_ca,'LYSO_CA','mmolCa/m2/d', 'Incorrect benthic output CaCO3', missing_value=0.0_rk)
   end subroutine initialize

   subroutine do_bottom(self, _ARGUMENTS_DO_)
      class(type_medusa_benthic), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: ZSEDC,ZSEDN,ZSEDFE,ZSEDSI,ZSEDCA,ZDET,ZDTC,ZSIL,ZALK,ZOXY
      real(rk), parameter :: s_per_d = 86400.0_rk
      real(rk) :: f_benout_c,f_benout_n,f_benout_fe,f_benout_ca,f_benout_si,cal_ccd,depth

      _BOTTOM_LOOP_BEGIN_

         _GET_BOTTOM_(self%id_ZSEDC,ZSEDC)
         _GET_BOTTOM_(self%id_ZSEDN,ZSEDN)
         _GET_BOTTOM_(self%id_ZSEDFE,ZSEDFE)
         _GET_BOTTOM_(self%id_ZSEDSI,ZSEDSI)
         _GET_BOTTOM_(self%id_ZSEDCA,ZSEDCA)
         _GET_(self%id_ZDET,ZDET)
         _GET_(self%id_ZDTC,ZDTC)
         _GET_(self%id_ZOXY,ZOXY)

         f_benout_c =  self%xsedc * ZSEDC
         f_benout_n =  self%xsedn * ZSEDN
         f_benout_fe = self%xsedfe * ZSEDFE
         f_benout_si = self%xsedsi * ZSEDSI
         f_benout_ca = self%xsedca * ZSEDCA

         ! organic components
         _ADD_BOTTOM_SOURCE_(self%id_ZSEDC, -f_benout_c)
         _ADD_BOTTOM_FLUX_(self%id_ZDIC,    +f_benout_c)

         _ADD_BOTTOM_SOURCE_(self%id_ZSEDN, -f_benout_n)
         _ADD_BOTTOM_FLUX_(self%id_ZDIN,    +f_benout_n)

         if (ZOXY .ge. self%xo2min) _ADD_BOTTOM_FLUX_(self%id_ZOXY, - self%xthetanit * f_benout_n - self%xthetarem * f_benout_c)

         _ADD_BOTTOM_SOURCE_(self%id_ZSEDFE, -f_benout_fe)
         _ADD_BOTTOM_FLUX_(self%id_ZFER,     +f_benout_n * self%xrfn)

         ! inorganic components
         _ADD_BOTTOM_SOURCE_(self%id_ZSEDSI, -f_benout_si)
         _ADD_BOTTOM_FLUX_(self%id_ZSIL,     +f_benout_si)

         _ADD_BOTTOM_SOURCE_(self%id_ZSEDCA, -f_benout_ca)
         _ADD_BOTTOM_FLUX_(self%id_ZALK,     + 2._rk * f_benout_ca)
         _ADD_BOTTOM_FLUX_(self%id_ZDIC,     +f_benout_ca)

         _SET_BOTTOM_DIAGNOSTIC_(self%id_f_benout_c, f_benout_c  * s_per_d)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_f_benout_n, f_benout_n  * s_per_d)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_f_benout_fe,f_benout_fe * s_per_d)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_f_benout_si,f_benout_si * s_per_d)
         _SET_BOTTOM_DIAGNOSTIC_(self%id_f_benout_ca,f_benout_ca * s_per_d)

         ! account for CaCO3 that dissolves when it shouldn't
         _GET_HORIZONTAL_(self%id_CAL_CCD,cal_ccd)
         _GET_(self%id_depth,depth)
         if (depth .le. cal_ccd) then   
            _SET_BOTTOM_DIAGNOSTIC_(self%id_f_lyso_ca,f_benout_ca * s_per_d)
         else
            _SET_BOTTOM_DIAGNOSTIC_(self%id_f_lyso_ca, 0._rk)
         end if

      _BOTTOM_LOOP_END_

   end subroutine do_bottom

end module medusa_benthic
