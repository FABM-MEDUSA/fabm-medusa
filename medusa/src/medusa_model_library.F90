module medusa_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use medusa_pelagic
   use medusa_oxygen

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: medusa_model_factory

contains

   subroutine create(self,name,model)
      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('medusa_pelagic');                          allocate(type_medusa_pelagic::model)
         case ('medusa_oxygen');                           allocate(type_medusa_oxygen::model)
         ! Add new models here
      end select
   end subroutine create

end module
