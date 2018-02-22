for d in $@
do
    includes=$(grep "#include" $d -Hr --include=*xx --include=*pp)
    for f in  boundary_factory.hxx boundary_op.hxx boundary_region.hxx \
                               boundary_standard.hxx boutcomm.hxx \
                               boutexception.hxx boutmain.hxx bout_types.hxx \
                               comm_group.hxx cyclic_reduction.hxx \
                               datafile.hxx dataformat.hxx dcomplex.hxx \
                               derivs.hxx difops.hxx fft.hxx field2d.hxx \
                               field3d.hxx field_data.hxx field_factory.hxx \
                               field.hxx fieldperp.hxx globals.hxx \
                               gyro_average.hxx initialprofiles.hxx \
                               interpolation_factory.hxx interpolation.hxx \
                               invert_laplace.hxx invert_parderiv.hxx \
                               lapack_routines.hxx mask.hxx msg_stack.hxx \
                               multiostream.hxx options.hxx \
                               optionsreader.hxx output.hxx \
                               parallel_boundary_op.hxx \
                               parallel_boundary_region.hxx smoothing.hxx \
                               sourcex.hxx stencils.hxx unused.hxx utils.hxx \
                               vecops.hxx vector2d.hxx vector3d.hxx where.hxx 
    do
        matched=$(echo "$includes"|grep [^/]$f)
        todo=$(echo "$matched"|cut -d: -f1)
        for doit in $todo
        do
            sed -i "s/\([<\"]\)$f/\1bout\/$f/" $doit
        done
    done
done
