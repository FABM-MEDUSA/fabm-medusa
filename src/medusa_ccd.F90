#include "fabm_driver.h"
module medusa_ccd

!************************************************************
! calculate compensation depth for calcite and aragonite
! providing corresponding saturation states as input
!************************************************************

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_medusa_ccd
      type (type_dependency_id)                     :: id_om1,id_dz
      type (type_horizontal_diagnostic_variable_id) :: id_CCD
   contains
      procedure :: initialize
      procedure :: do_column
   end type

contains

   subroutine initialize(self, configunit)
      class (type_medusa_ccd), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      call self%register_implemented_routines((/source_do_column/))
      call self%register_diagnostic_variable(self%id_CCD,'CCD','m','CCD depth',source=source_do_column)
      call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
      call self%register_dependency(self%id_om1,'OMEGA','-','Omega calcite 3D')
   end subroutine

   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class(type_medusa_ccd), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

       real(rk) :: depth,depthb,dz,f_om,f_om_a
       real(rk) :: fq0,fq1,fq2,fq3,fq4, f2_ccd_cal

       depth = 0._rk
       depthb = 0._rk

      _DOWNWARD_LOOP_BEGIN_

         _GET_(self%id_dz, dz)
         _GET_(self%id_om1, f_om)

         depth = depth + 0.5_rk * dz

         if (f_om .lt. 1._rk) then
            if (depthb == 0._rk) then
               ! In surface layer (use centre depth)
               f2_ccd_cal = depth
            else
               ! Below surface layer (linear interpolation to find depth at which f_om==1)
               fq0 = f_om_a - f_om
               fq1 = f_om_a - 1._rk
               fq2 = fq1 / (fq0 + tiny(fq0))
               fq3 = depth-depthb
               fq4 = fq2 * fq3
               f2_ccd_cal = depthb + fq4
            endif
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_CCD, f2_ccd_cal)
            return
         endif

         depthb = depth
         f_om_a = f_om

         depth = depth + 0.5_rk * dz

      _DOWNWARD_LOOP_END_

      ! reached seafloor and still no dissolution; set to seafloor depth (W-point)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_CCD, depth)

   end subroutine do_column

end module
