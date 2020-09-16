function mpi_update_fields(grid, new_field, old_field, 
                           new_proc_to_blocks_mapping, old_proc_to_blocks_mapping,
                           new_block_to_proc_mapping,  old_block_to_proc_mapping )

    i_proc = MPI.Comm_rank(MPI.COMM_WORLD) + 1
    new_my_blocks = new_proc_to_blocks_mapping[i_proc]
    old_my_blocks = old_proc_to_blocks_mapping[i_proc]
    println("new_proc_to_blocks_mapping: $new_proc_to_blocks_mapping")
    println("old_proc_to_blocks_mapping: $old_proc_to_blocks_mapping")
    println("new_block_to_proc_mapping: $new_block_to_proc_mapping")
    println("old_block_to_proc_mapping: $old_block_to_proc_mapping")

    n_points = (2*grid.ngx+grid.nx) * (2*grid.ngy+grid.ny) * (2*grid.ngz+grid.nz)

    requests = Array{MPI.Request,1}()

    # Send blocks that I had (in old_my_blocks) but I do not have any longer (new_proc != i_proc)
    for (ib, my_ib) in old_my_blocks
        if ib > 0
            new_proc = new_block_to_proc_mapping[ib][1]
            if new_proc != i_proc
                println("mpi case 1 : $ib")
                tag = ib
                request = MPI.Isend(old_field[my_ib].val, new_proc-1, tag, MPI.COMM_WORLD)
                push!(requests, request)
            end
        end
    end

    for (ib, my_ib) in new_my_blocks
        if ib > 0
            old_proc = old_block_to_proc_mapping[ib][1]
            # Copy blocks that I have (in new_my_blocks) and I had, too (old_proc == i_proc)
            if old_proc == i_proc
                println("mpi case 2 : $ib")
                old_pos = old_block_to_proc_mapping[ib][2]
                println("my_ib, old_pos, i_proc: $my_ib, $old_pos, $i_proc")
                new_field[my_ib].val[:,:,:] = old_field[old_pos].val[:,:,:]  # [:,:,:] are required to deep copy
                # Receive blocks that I have (in new_my_blocks) but I did not have (old_proc != i_proc)
            else
                println("mpi case 3 : $ib")
                tag = ib
                request = MPI.Irecv!(new_field[my_ib].val, old_proc-1, tag, MPI.COMM_WORLD)
                push!(requests, request)
            end
        end
    end

    stats = MPI.Waitall!(requests)

end
