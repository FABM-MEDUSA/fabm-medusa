#include "fabm_driver.h"
module medusa_ccd

!************************************************************
! calculate compensation depth for calcite and aragonite
! providing corresponding saturation states as input
!************************************************************

   use fabm_types
   use fabm_standard_variables

   implicit none

   private

   type,extends(type_base_model),public :: type_medusa_ccd

      type (type_dependency_id)            :: id_depth,id_om1,id_dz
      type (type_horizontal_diagnostic_variable_id) :: id_CCD

   contains
     procedure :: initialize
     procedure :: get_light => do_ccd

   end type

contains
    
    subroutine initialize(self,configunit)

      class (type_medusa_ccd), intent(inout), target :: self
      integer,                      intent(in)            :: configunit

     call self%register_diagnostic_variable(self%id_CCD,'CCD','m','CCD depth',source=source_do_column)
     call self%register_dependency(self%id_dz, standard_variables%cell_thickness)
     call self%register_dependency(self%id_depth, standard_variables%depth)
     call self%register_dependency(self%id_om1,'OMEGA','-','Omega calcite 3D')

    end subroutine

   subroutine do_ccd(self,_ARGUMENTS_VERTICAL_)
   class(type_medusa_ccd),intent(in) :: self
   
   _DECLARE_ARGUMENTS_VERTICAL_

    real(rk) :: depth,depthb,dz,f_om,f_om_a
    real(rk) :: fq0,fq1,fq2,fq3,fq4, f2_ccd_cal
    integer :: i2_om,check

    i2_om = 0
    check = 0 !0 for surface layer

   _VERTICAL_LOOP_BEGIN_

     _GET_(self%id_dz,dz)
     _GET_(self%id_om1, f_om)
     _GET_(self%id_depth,depth)

    if ((i2_om == 0).and.(f_om .lt. 1._rk)) then
       if (check == 0) then 
           f2_ccd_cal = depth
       else
          fq0 = f_om_a - f_om
          fq1 = f_om_a - 1._rk
          fq2 = fq1 / (fq0 + tiny(fq0))
          fq3 = depth-depthb
          fq4 = fq2 * fq3
          f2_ccd_cal = depthb + fq4
       endif
    i2_om = 1
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_CCD, f2_ccd_cal)
    endif

    check = 1
    depthb = depth
    f_om_a = f_om
   _VERTICAL_LOOP_END_

   if (i2_om == 0) then ! reached seafloor and still no dissolution; set to seafloor (W-point)
    _SET_HORIZONTAL_DIAGNOSTIC_(self%id_CCD, depth + 0.5_rk * dz)
    endif

   end subroutine do_ccd

end module
