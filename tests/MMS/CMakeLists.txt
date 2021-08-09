add_subdirectory(advection)
add_subdirectory(diffusion)
add_subdirectory(diffusion2)
add_subdirectory(hw)
add_subdirectory(laplace)
add_subdirectory(spatial/advection)
add_subdirectory(spatial/d2dx2)
add_subdirectory(spatial/d2dz2)
add_subdirectory(spatial/diffusion)
add_subdirectory(spatial/fci)
add_subdirectory(time)
add_subdirectory(time-petsc)
add_subdirectory(wave-1d)
add_subdirectory(wave-1d-y)

########################################
# The following require boutcore:

# add_subdirectory(bracket)
# add_subdirectory(derivatives3)
# add_subdirectory(shiftedmetricinterp)
# add_subdirectory(upwinding3)

########################################
# The following are marked as broken

# add_subdirectory(elm-pb)
# add_subdirectory(fieldalign)
# add_subdirectory(GBS)
# add_subdirectory(tokamak)
