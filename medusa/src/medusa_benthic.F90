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
      type (type_state_variable_id)        ::
      type (type_dependency_id)            ::
      type (type_diagnostic_variable_id)   ::
      type (type_state_variable_id)   ::
      ! Parameters
      logical ::
      real(rk) :: 


   contains

      procedure :: initialize
      procedure :: do_bottom

  end type type_medusa_benthic

contains

   subroutine initialize(self,configunit)

   class(type_medusa_benthic),intent(inout),target :: self
   integer,               intent(in)           :: configunit
   real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_)

   class(type_medusa_benthic), INTENT(IN) :: self
  _DECLARE_ARGUMENTS_DO_BOTTOM_

! !LOCAL VARIABLES:

    real(rk) ::
    real(rk), parameter :: d_per_s = 1.0_rk/86400.0_rk

    _HORIZONTAL_LOOP_BEGIN_

   _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

  end module medusa_benthic
