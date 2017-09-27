add_library(fabm_models_medusa OBJECT
            medusa_model_library.F90
            medusa_pelagic.F90
            medusa_fast_detritus.F90
            medusa_oxygen.F90
            medusa_carbonate.F90
           )

add_dependencies(fabm_models_medusa fabm_base)
