function refine_boundary_minx(grid, f::Array{Float64,3}, iy::Int64, iz::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(ngx, ny, nz)
    # loop over coarse grid
    for ka=iz*div(nz,2)+1:(iz+1)*div(nz,2), ja=iy*div(ny,2)+1:(iy+1)*div(ny,2), ia=1:div(ngx,2)
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx-ngx; j_ref = 2*ja-1+ngy-iy*ny-ngy; k_ref = 2*ka-1+ngz-iz*nz-ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function refine_boundary_maxx(grid, f::Array{Float64,3}, iy::Int64, iz::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(ngx, ny, nz)
    # loop over coarse grid
    for ka=iz*div(nz,2)+1:(iz+1)*div(nz,2), ja=iy*div(ny,2)+1:(iy+1)*div(ny,2), ia=nx-div(ngx,2)+1:nx
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1-2*nx+ngx; j_ref = 2*ja-1+ngy-iy*ny-ngy; k_ref = 2*ka-1+ngz-iz*nz-ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function derefine_boundary_minx(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(ngx, div(ny,2), div(nz,2))
    # loop over coarse grid
    for ka=1:div(nz,2), ja=1:div(ny,2), ia=1:ngx
        i=ia;           j=ja;           k=ka  # no ng* added here, boundaries are without corners
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function derefine_boundary_maxx(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(ngx, div(ny,2), div(nz,2))
    # loop over coarse grid
    for ka=1:div(nz,2), ja=1:div(ny,2), ia=div(nx,2)-ngx+1:div(nx,2)
        i=ia-div(nx,2)+ngx;           j=ja;           k=ka
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function refine_boundary_miny(grid, f::Array{Float64,3}, ix::Int64, iz::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(nx, ngy, nz)
    # loop over coarse grid
    for ka=iz*div(nz,2)+1:(iz+1)*div(nz,2), ja=1:div(ngy,2), ia=ix*div(nx,2)+1:(ix+1)*div(nx,2)
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx-ix*nx-ngx; j_ref = 2*ja-1+ngy-ngy; k_ref = 2*ka-1+ngz-iz*nz-ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function refine_boundary_maxy(grid, f::Array{Float64,3}, ix::Int64, iz::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(nx, ngy, nz)
    # loop over coarse grid
    for ka=iz*div(nz,2)+1:(iz+1)*div(nz,2), ja=ny-div(ngy,2)+1:ny, ia=ix*div(nx,2)+1:(ix+1)*div(nx,2)
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx-ix*nx-ngx; j_ref = 2*ja-1-2*ny+ngy; k_ref = 2*ka-1+ngz-iz*nz-ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function derefine_boundary_miny(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(div(nx,2), ngy, div(nz,2))
    # loop over coarse grid
    for ka=1:div(nz,2), ja=1:ngy, ia=1:div(nx,2)
        i=ia;           j=ja;           k=ka  # no ng* added here, boundaries are without corners
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function derefine_boundary_maxy(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(div(nx,2), ngy, div(nz,2))
    # loop over coarse grid
    for ka=1:div(nz,2), ja=div(ny,2)-ngy+1:div(ny,2), ia=1:div(nx,2)
        i=ia ;  j=ja-div(ny,2)+ngy;  k=ka
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function refine_boundary_minz(grid, f::Array{Float64,3}, ix::Int64, iy::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(nx, ny, ngz)
    # loop over coarse grid
    for ka=1:div(ngz,2), ja=iy*div(ny,2)+1:(iy+1)*div(ny,2), ia=ix*div(nx,2)+1:(ix+1)*div(nx,2)
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx-ix*nx-ngx; j_ref = 2*ja-1+ngy-iy*ny-ngy; k_ref = 2*ka-1+ngz-ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function refine_boundary_maxz(grid, f::Array{Float64,3}, ix::Int64, iy::Int64)
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_ref = zeros(nx, ny, ngz)
    # loop over coarse grid
    for ka=nz-div(ngz,2)+1:nz, ja=iy*div(ny,2)+1:(iy+1)*div(ny,2), ia=ix*div(nx,2)+1:(ix+1)*div(nx,2)
        i=ia+ngx;           j=ja+ngy;           k=ka+ngz
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx-ix*nx-ngx; j_ref = 2*ja-1+ngy-iy*ny-ngy; k_ref = 2*ka-1-2*nz+ngz; 
        # refine from 3*3*3 coarse points to 2*2*2 fine points surrounding the coarse central point
        #f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = refine_points(f[i-1:i+1, j-1:j+1, k-1:k+1])
        f_ref[i_ref:i_ref+1, j_ref:j_ref+1, k_ref:k_ref+1] = f[i,j,k]
    end
    return f_ref
end

function derefine_boundary_minz(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(div(nx,2), div(ny,2), ngz)
    # loop over coarse grid
    for ka=1:ngz, ja=1:div(ny,2), ia=1:div(nx,2)
        i=ia;           j=ja;           k=ka  # no ng* added here, boundaries are without corners
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function derefine_boundary_maxz(grid, f::Array{Float64,3})
    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz
    f_deref = zeros(div(nx,2), div(ny,2), ngz)
    # loop over coarse grid
    for ka=div(nz,2)-ngz+1:div(nz,2), ja=1:div(ny,2), ia=1:div(nx,2)
        i=ia ; j=ja ;  k=ka-div(nz,2)+ngz;
        # get the corresponding fine grid point "below"
        # indexes are translated to have the start index equal to 1
        i_ref = 2*ia-1+ngx; j_ref = 2*ja-1+ngy; k_ref = 2*ka-1+ngz; 
        # derefine from 3*3*3 fine points (i_ref:i_ref+2 could be changed to i_ref-1:i_ref+1)
        #f_deref[i, j, k] = derefine_point(f[i_ref:i_ref+2, j_ref:j_ref+2, k_ref:k_ref+2])
        f_deref[i, j, k] = 0.125*(f[i_ref,j_ref,k_ref]     + f[i_ref+1,j_ref,k_ref]   + 
                                  f[i_ref,j_ref+1,k_ref]   + f[i_ref+1,j_ref+1,k_ref] +
                                  f[i_ref,j_ref,k_ref+1]   + f[i_ref+1,j_ref,k_ref+1] +
                                  f[i_ref,j_ref+1,k_ref+1] + f[i_ref+1,j_ref+1,k_ref+1])
    end
    return f_deref
end

function mpi_update_ghost(grid, field_name, proc_to_blocks_mapping, block_to_proc_mapping)

    i_proc = MPI.Comm_rank(MPI.COMM_WORLD) + 1
    n_procs = MPI.Comm_size(MPI.COMM_WORLD)  
    n_blocks = size(grid.blocks,1)
    my_blocks = proc_to_blocks_mapping[i_proc]
    println("GHOST proc_to_blocks_mapping: $proc_to_blocks_mapping")
    println("GHOST block_to_proc_mapping:  $block_to_proc_mapping")

    n_points = (2*grid.ngx+grid.nx) * (2*grid.ngy+grid.ny) * (2*grid.ngz+grid.nz)

    requests = Array{MPI.Request,1}()

    nx = grid.nx; ny = grid.ny; nz = grid.nz; ngx = grid.ngx; ngy = grid.ngy; ngz = grid.ngz

    # Since interpolation (refine/derefine) uses 1-level ghost nodes I here update ghost 
    # edges and corners extrapolating values. It does not make sense much because
    # here I am udpating ghosts while here I use to values to extrapolate the values
    # theirselves. Probably it would be better to use ghost interpolation which do
    # not rely on ghost values theirsevles!
    #NOT NEEDED WITH CURRENT INTERPOLATIONS for (ib, my_ib) in proc_to_blocks_mapping[i_proc]
    #NOT NEEDED WITH CURRENT INTERPOLATIONS     println("fixing block edge and corners $ib")
    #NOT NEEDED WITH CURRENT INTERPOLATIONS     fix_edge_and_corners(grid, field_name[my_ib].val)
    #NOT NEEDED WITH CURRENT INTERPOLATIONS end

    #-----------------------------------------------------------------------------
    # I need buffers for different reasons:
    #-----------------------------------------------------------------------------
    # 1) to match the send/recv and for performance reasons it is far better to
    #    have proc i to send/recv proc j only one (or a few buffers) not each
    #    block face because a processor can have many many blocks
    # 2) recv cannot work if passing array sections, therefore buffers on 
    #    receiving are mandatory
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # [1] Get the list of procs I need to talk to exchange boundaries
    #-----------------------------------------------------------------------------
    procs_adj_blocks_minx_ref0 = Array{Int64,1}()
    procs_adj_blocks_maxx_ref0 = Array{Int64,1}()
    procs_adj_blocks_minx_refp1 = Array{Int64,1}()
    procs_adj_blocks_maxx_refp1 = Array{Int64,1}()
    procs_adj_blocks_minx_refm1 = Array{Int64,1}()
    procs_adj_blocks_maxx_refm1 = Array{Int64,1}()

    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        bc = block.boundaries[1]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
    #        println("KIT i_proc: $i_proc - ib: $ib - adj_ib: $adj_ib - owner_adj_ib: $owner_adj_ib")
            if owner_adj_ib != i_proc
                if !in(owner_adj_ib, procs_adj_blocks_minx_ref0)
                    push!(procs_adj_blocks_minx_ref0, owner_adj_ib)
                end
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                if !in(owner_adj_ib, procs_adj_blocks_minx_refp1)
                    push!(procs_adj_blocks_minx_refp1, owner_adj_ib)
                end
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                if owner_adj_ib != i_proc
                    if !in(owner_adj_ib, procs_adj_blocks_minx_refm1)
                        push!(procs_adj_blocks_minx_refm1, owner_adj_ib)
                    end
                end
            end
        end

        # xmax boundary
        bc = block.boundaries[2]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
    #        println("KIT i_proc: $i_proc - ib: $ib - adj_ib: $adj_ib - owner_adj_ib: $owner_adj_ib")
            if owner_adj_ib != i_proc
                if !in(owner_adj_ib, procs_adj_blocks_maxx_ref0)
                    push!(procs_adj_blocks_maxx_ref0, owner_adj_ib)
                end
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                if !in(owner_adj_ib, procs_adj_blocks_maxx_refp1)
                    push!(procs_adj_blocks_maxx_refp1, owner_adj_ib)
                end
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                if owner_adj_ib != i_proc
                    if !in(owner_adj_ib, procs_adj_blocks_maxx_refm1)
                        push!(procs_adj_blocks_maxx_refm1, owner_adj_ib)
                    end
                end
            end
        end
    end

    n_procs_adj_blocks_minx_ref0 = size(procs_adj_blocks_minx_ref0, 1)
    n_procs_adj_blocks_maxx_ref0 = size(procs_adj_blocks_maxx_ref0, 1)
    n_procs_adj_blocks_minx_refp1 = size(procs_adj_blocks_minx_refp1, 1)
    n_procs_adj_blocks_maxx_refp1 = size(procs_adj_blocks_maxx_refp1, 1)
    n_procs_adj_blocks_minx_refm1 = size(procs_adj_blocks_minx_refm1, 1)
    n_procs_adj_blocks_maxx_refm1 = size(procs_adj_blocks_maxx_refm1, 1)

    println("GHOST1 minx i_proc, procs_adj_blocks_minx_ref0:  $i_proc => $(procs_adj_blocks_minx_ref0)")
    println("GHOST1 maxx i_proc, procs_adj_blocks_maxx_ref0:  $i_proc => $(procs_adj_blocks_maxx_ref0)")
    println("GHOST1 minx i_proc, procs_adj_blocks_minx_refp1: $i_proc => $(procs_adj_blocks_minx_refp1)")
    println("GHOST1 maxx i_proc, procs_adj_blocks_maxx_refp1: $i_proc => $(procs_adj_blocks_maxx_refp1)")
    println("GHOST1 minx i_proc, procs_adj_blocks_minx_refm1: $i_proc => $(procs_adj_blocks_minx_refm1)")
    println("GHOST1 maxx i_proc, procs_adj_blocks_maxx_refm1: $i_proc => $(procs_adj_blocks_maxx_refm1)")

    #-----------------------------------------------------------------------------
    # [2] Prepare list of blocks to send/recv for each proc
    #-----------------------------------------------------------------------------
    blocks_per_proc_send_minx_ref0 = Array{Any,1}()
    blocks_per_proc_recv_minx_ref0 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_minx_ref0
        push!(blocks_per_proc_send_minx_ref0, Array{Int64,1}())
        push!(blocks_per_proc_recv_minx_ref0, Array{Int64,1}())
    end
    blocks_per_proc_send_maxx_ref0 = Array{Any,1}()
    blocks_per_proc_recv_maxx_ref0 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_maxx_ref0
        push!(blocks_per_proc_send_maxx_ref0, Array{Int64,1}())
        push!(blocks_per_proc_recv_maxx_ref0, Array{Int64,1}())
    end

    blocks_per_proc_send_minx_refp1 = Array{Any,1}()
    blocks_per_proc_recv_minx_refp1 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_minx_refp1
        push!(blocks_per_proc_send_minx_refp1, Array{Int64,1}())
        push!(blocks_per_proc_recv_minx_refp1, Array{Int64,1}())
    end
    blocks_per_proc_send_maxx_refp1 = Array{Any,1}()
    blocks_per_proc_recv_maxx_refp1 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_maxx_refp1
        push!(blocks_per_proc_send_maxx_refp1, Array{Int64,1}())
        push!(blocks_per_proc_recv_maxx_refp1, Array{Int64,1}())
    end

    blocks_per_proc_send_minx_refm1 = Array{Any,1}()
    blocks_per_proc_recv_minx_refm1 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_minx_refm1
        push!(blocks_per_proc_send_minx_refm1, Array{Any,1}()) # tuple with block and 1 to 4 pos
        push!(blocks_per_proc_recv_minx_refm1, Array{Any,1}()) # tuple with block and 1 to 4 pos
    end
    blocks_per_proc_send_maxx_refm1 = Array{Any,1}()
    blocks_per_proc_recv_maxx_refm1 = Array{Any,1}()
    for i = 1:n_procs_adj_blocks_maxx_refm1
        push!(blocks_per_proc_send_maxx_refm1, Array{Any,1}()) # tuple with block and 1 to 4 pos
        push!(blocks_per_proc_recv_maxx_refm1, Array{Any,1}()) # tuple with block and 1 to 4 pos
    end
    println("GHOST2 done")

    #-----------------------------------------------------------------------------
    # [3] Build list of blocks to send for each proc
    #-----------------------------------------------------------------------------
    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        bc = block.boundaries[1]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_ref0, owner_adj_ib)
                push!(blocks_per_proc_send_minx_ref0[pos_adj_owner_proc], ib)
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_refp1, owner_adj_ib)
                push!(blocks_per_proc_send_minx_refp1[pos_adj_owner_proc], ib)
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                if owner_adj_ib != i_proc
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_refm1, owner_adj_ib)
                    push!(blocks_per_proc_send_minx_refm1[pos_adj_owner_proc], (ib, i_adj_ib))
                end
            end
        end

        # xmax boundary
        bc = block.boundaries[2]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_ref0, owner_adj_ib)
                push!(blocks_per_proc_send_maxx_ref0[pos_adj_owner_proc], ib)
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib != i_proc
                pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_refp1, owner_adj_ib)
                push!(blocks_per_proc_send_maxx_refp1[pos_adj_owner_proc], ib)
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                if owner_adj_ib != i_proc
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_refm1, owner_adj_ib)
                    push!(blocks_per_proc_send_maxx_refm1[pos_adj_owner_proc], (ib, i_adj_ib))
                end
            end
        end
    end

    println("GHOST3 i_proc, blocks_per_proc_send_minx_ref0:  $i_proc => $(blocks_per_proc_send_minx_ref0)")
    println("GHOST3 i_proc, blocks_per_proc_send_maxx_ref0:  $i_proc => $(blocks_per_proc_send_maxx_ref0)")
    println("GHOST3 i_proc, blocks_per_proc_send_minx_refp1: $i_proc => $(blocks_per_proc_send_minx_refp1)")
    println("GHOST3 i_proc, blocks_per_proc_send_maxx_refp1: $i_proc => $(blocks_per_proc_send_maxx_refp1)")
    println("GHOST3 i_proc, blocks_per_proc_send_minx_refm1: $i_proc => $(blocks_per_proc_send_minx_refm1)")
    println("GHOST3 i_proc, blocks_per_proc_send_maxx_refm1: $i_proc => $(blocks_per_proc_send_maxx_refm1)")

    #-----------------------------------------------------------------------------
    # [4] Getting the blocks to recv from 
    #-----------------------------------------------------------------------------
    for proc_recv_from in procs_adj_blocks_minx_ref0
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_ref0=$procs_adj_blocks_minx_ref0")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmax boundary
            bc = block.boundaries[2]
            if bc.ref_change == 0
                adj_ib = bc.adj_ibs[1]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                # println("KKK i_proc=$i_proc - proc_recv_from=$proc_recv_from - ib=$ib - my_ib=$my_ib - owner_adj_ib=$owner_adj_ib")
                if owner_adj_ib == i_proc
                    # TODO basta mettere enumerate nel loop fuori
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_ref0, proc_recv_from)
                    push!(blocks_per_proc_recv_minx_ref0[pos_adj_owner_proc], adj_ib)
                end
            end
        end
    end

    for proc_recv_from in procs_adj_blocks_maxx_ref0
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_ref0=$procs_adj_blocks_minx_ref0")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmin boundary
            bc = block.boundaries[1]
            if bc.ref_change == 0
                adj_ib = bc.adj_ibs[1]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                if owner_adj_ib == i_proc
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_ref0, proc_recv_from)
                    push!(blocks_per_proc_recv_maxx_ref0[pos_adj_owner_proc], adj_ib)
                end
            end
        end
    end

    for proc_recv_from in procs_adj_blocks_minx_refp1
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_refp1=$procs_adj_blocks_minx_refp1")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmax boundary
            bc = block.boundaries[2]
            if bc.ref_change == -1
                for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                    owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                    # println("KKK i_proc=$i_proc - proc_recv_from=$proc_recv_from - ib=$ib - my_ib=$my_ib - owner_adj_ib=$owner_adj_ib")
                    if owner_adj_ib == i_proc
                        # TODO basta mettere enumerate nel loop fuori
                        pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_refp1, proc_recv_from)
                        push!(blocks_per_proc_recv_minx_refp1[pos_adj_owner_proc], adj_ib)
                    end
                end
            end
        end
    end

    for proc_recv_from in procs_adj_blocks_maxx_refp1
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_refp1=$procs_adj_blocks_minx_refp1")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmin boundary
            bc = block.boundaries[1]
            if bc.ref_change == -1
                for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                    owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                    if owner_adj_ib == i_proc
                        pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_refp1, proc_recv_from)
                        push!(blocks_per_proc_recv_maxx_refp1[pos_adj_owner_proc], adj_ib)
                    end
                end
            end
        end
    end

    for proc_recv_from in procs_adj_blocks_minx_refm1
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_refm1=$procs_adj_blocks_minx_refm1")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmax boundary
            bc = block.boundaries[2]
            if bc.ref_change == 1
                adj_ib = bc.adj_ibs[1]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                min_y = block.min_coord_max_ref_y
                min_z = block.min_coord_max_ref_z
                min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
                min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
                if min_z == min_z_adj
                    if min_y == min_y_adj
                        pos_1_4 = 1
                    else
                        pos_1_4 = 2
                    end
                else
                    if min_y == min_y_adj
                        pos_1_4 = 3
                    else
                        pos_1_4 = 4
                    end
                end
                # println("KKK i_proc=$i_proc - proc_recv_from=$proc_recv_from - ib=$ib - my_ib=$my_ib - owner_adj_ib=$owner_adj_ib")
                if owner_adj_ib == i_proc
                    # TODO basta mettere enumerate nel loop fuori
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_refm1, proc_recv_from)
                    push!(blocks_per_proc_recv_minx_refm1[pos_adj_owner_proc], (adj_ib, pos_1_4))
                end
            end
        end
    end

    for proc_recv_from in procs_adj_blocks_maxx_refm1
        println("KAI i_proc=$i_proc - proc_recv_from=$proc_recv_from - procs_adj_blocks_minx_refm1=$procs_adj_blocks_minx_refm1")
        for (ib, my_ib) in proc_to_blocks_mapping[proc_recv_from]
            block = grid.blocks[ib]

            # xmin boundary
            bc = block.boundaries[1]
            if bc.ref_change == 1
                adj_ib = bc.adj_ibs[1]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                min_y = block.min_coord_max_ref_y
                min_z = block.min_coord_max_ref_z
                min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
                min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
                if min_z == min_z_adj
                    if min_y == min_y_adj
                        pos_1_4 = 1
                    else
                        pos_1_4 = 2
                    end
                else
                    if min_y == min_y_adj
                        pos_1_4 = 3
                    else
                        pos_1_4 = 4
                    end
                end
                if owner_adj_ib == i_proc
                    pos_adj_owner_proc = findfirst(procs_adj_blocks_maxx_refm1, proc_recv_from)
                    push!(blocks_per_proc_recv_maxx_refm1[pos_adj_owner_proc], (adj_ib, pos_1_4))
                end
            end
        end
    end

    println("GHOST4 i_proc, blocks_per_proc_recv_minx_ref0: $i_proc => $(blocks_per_proc_recv_minx_ref0)")
    println("GHOST4 i_proc, blocks_per_proc_recv_maxx_ref0: $i_proc => $(blocks_per_proc_recv_maxx_ref0)")
    println("GHOST4 i_proc, blocks_per_proc_recv_minx_refp1: $i_proc => $(blocks_per_proc_recv_minx_refp1)")
    println("GHOST4 i_proc, blocks_per_proc_recv_maxx_refp1: $i_proc => $(blocks_per_proc_recv_maxx_refp1)")
    println("GHOST4 i_proc, blocks_per_proc_recv_minx_refm1: $i_proc => $(blocks_per_proc_recv_minx_refm1)")
    println("GHOST4 i_proc, blocks_per_proc_recv_maxx_refm1: $i_proc => $(blocks_per_proc_recv_maxx_refm1)")

    #-----------------------------------------------------------------------------
    # [5] Preparing the buffers to send/recv
    #-----------------------------------------------------------------------------
    buffers_send_minx_ref0 = Array{Any,1}()
    buffers_recv_minx_ref0 = Array{Any,1}()
    buffers_send_maxx_ref0 = Array{Any,1}()
    buffers_recv_maxx_ref0 = Array{Any,1}()
    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        n_blocks_per_proc = size(blocks_per_proc_send_minx_ref0[pos_proc],1)
        push!( buffers_send_minx_ref0, zeros(ngx, ny, nz, n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_minx_ref0[pos_proc],1)
        push!( buffers_recv_minx_ref0, zeros(ngx, ny, nz, n_blocks_per_proc) )
    end
    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_ref0)
        n_blocks_per_proc = size(blocks_per_proc_send_maxx_ref0[pos_proc],1)
        push!( buffers_send_maxx_ref0, zeros(ngx, ny, nz, n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_maxx_ref0[pos_proc],1)
        push!( buffers_recv_maxx_ref0, zeros(ngx, ny, nz, n_blocks_per_proc) )
    end

    buffers_send_minx_refp1 = Array{Any,1}()
    buffers_recv_minx_refp1 = Array{Any,1}()
    buffers_send_maxx_refp1 = Array{Any,1}()
    buffers_recv_maxx_refp1 = Array{Any,1}()
    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refp1)
        n_blocks_per_proc = size(blocks_per_proc_send_minx_refp1[pos_proc],1)
        push!( buffers_send_minx_refp1, zeros(ngx, div(ny,2), div(nz,2), n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_minx_refp1[pos_proc],1)
        push!( buffers_recv_minx_refp1, zeros(ngx, ny, nz, n_blocks_per_proc) )
    end
    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refp1)
        n_blocks_per_proc = size(blocks_per_proc_send_maxx_refp1[pos_proc],1)
        push!( buffers_send_maxx_refp1, zeros(ngx, div(ny,2), div(nz,2), n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_maxx_refp1[pos_proc],1)
        push!( buffers_recv_maxx_refp1, zeros(ngx, ny, nz, n_blocks_per_proc) )
    end

    buffers_send_minx_refm1 = Array{Any,1}()
    buffers_recv_minx_refm1 = Array{Any,1}()
    buffers_send_maxx_refm1 = Array{Any,1}()
    buffers_recv_maxx_refm1 = Array{Any,1}()
    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refm1)
        n_blocks_per_proc = size(blocks_per_proc_send_minx_refm1[pos_proc],1)
        push!( buffers_send_minx_refm1, zeros(ngx, ny, nz, n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_minx_refm1[pos_proc],1)
        push!( buffers_recv_minx_refm1, zeros(ngx, div(ny,2), div(nz,2), n_blocks_per_proc) )
    end
    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refm1)
        n_blocks_per_proc = size(blocks_per_proc_send_maxx_refm1[pos_proc],1)
        push!( buffers_send_maxx_refm1, zeros(ngx, ny, nz, n_blocks_per_proc) )
        n_blocks_per_proc = size(blocks_per_proc_recv_maxx_refm1[pos_proc],1)
        push!( buffers_recv_maxx_refm1, zeros(ngx, div(ny,2), div(nz,2), n_blocks_per_proc) )
    end

    println("GHOST5 done")

    ##################################################################################

    #-----------------------------------------------------------------------------
    # [6] Filling buffers, starting SEND and RECV
    #-----------------------------------------------------------------------------
    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        for (pos_block, block) in enumerate(blocks_per_proc_send_minx_ref0[pos_proc])
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_minx_ref0[pos_proc][:, :, :, pos_block] = 
                field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
        end
        tag = 100 # send minx to maxx ref0
        request = MPI.Isend(buffers_send_minx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_ref0)
        tag = 100 # send minx to maxx ref0
        request = MPI.Irecv!(buffers_recv_maxx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refp1)
        for (pos_block, block) in enumerate(blocks_per_proc_send_minx_refp1[pos_proc])
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_minx_refp1[pos_proc][:, :, :, pos_block] = 
                derefine_boundary_minx(grid, field_name[my_ib].val)
        end
        tag = 110 # send minx to maxx refp1
        request = MPI.Isend(buffers_send_minx_refp1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refm1)
        tag = 110 # send minx to maxx refp1
        request = MPI.Irecv!(buffers_recv_maxx_refm1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refm1)
        for (pos_block, block_tuple) in enumerate(blocks_per_proc_send_minx_refm1[pos_proc])
            block = block_tuple[1]
            pos_1_4 = block_tuple[2]
            iy = [0, 1, 0, 1][pos_1_4]
            iz = [0, 0, 1, 1][pos_1_4]
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_minx_refm1[pos_proc][:, :, :, pos_block] = 
                refine_boundary_minx(grid, field_name[my_ib].val, iy, iz)
        end
        tag = 120 # send minx to maxx refm1
        request = MPI.Isend(buffers_send_minx_refm1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refp1)
        tag = 120 # send minx to maxx refm1
        request = MPI.Irecv!(buffers_recv_maxx_refp1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    #----------------------------------------------------------------------------------------

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_ref0)
        for (pos_block, block) in enumerate(blocks_per_proc_send_maxx_ref0[pos_proc])
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_maxx_ref0[pos_proc][:, :, :, pos_block] = 
                field_name[my_ib].val[nx+1:nx+ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
        end
        tag = 200 # send maxx to minx
        request = MPI.Isend(buffers_send_maxx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        tag = 200 # send maxx to minx
        request = MPI.Irecv!(buffers_recv_minx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refp1)
        for (pos_block, block) in enumerate(blocks_per_proc_send_maxx_refp1[pos_proc])
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_maxx_refp1[pos_proc][:, :, :, pos_block] = 
                derefine_boundary_maxx(grid, field_name[my_ib].val)
        end
        tag = 210 # send minx to maxx refp1
        request = MPI.Isend(buffers_send_maxx_refp1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refm1)
        tag = 210 # send minx to maxx refp1
        request = MPI.Irecv!(buffers_recv_minx_refm1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refm1)
        for (pos_block, block_tuple) in enumerate(blocks_per_proc_send_maxx_refm1[pos_proc])
            block = block_tuple[1]
            pos_1_4 = block_tuple[2]
            iy = [0, 1, 0, 1][pos_1_4]
            iz = [0, 0, 1, 1][pos_1_4]
            my_ib = block_to_proc_mapping[block][2]
            buffers_send_maxx_refm1[pos_proc][:, :, :, pos_block] = 
                refine_boundary_maxx(grid, field_name[my_ib].val, iy, iz)
        end
        tag = 220 # send minx to maxx refm1
        request = MPI.Isend(buffers_send_maxx_refm1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refp1)
        tag = 220 # send minx to maxx refm1
        request = MPI.Irecv!(buffers_recv_minx_refp1[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end
    println("GHOST6 done")

    #-----------------------------------------------------------------------------
    # [7] Send/recv minx/maxx from the same proc
    #-----------------------------------------------------------------------------
    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        bc = block.boundaries[1]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                    field_name[pos_adj_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            min_y = block.min_coord_max_ref_y
            min_z = block.min_coord_max_ref_z
            min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
            min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
            if owner_adj_ib == i_proc
                if min_y == min_y_adj
                    if min_z == min_z_adj
                        field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 0, 0)
                    else
                        field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 0, 1)
                    end
                else
                    if min_z == min_z_adj
                        field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 1, 0)
                    else
                        field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 1, 1)
                    end
                end
            end
        else # ref_change == -1
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                iy = [0, 1, 0, 1][i_adj_ib]
                iz = [0, 0, 1, 1][i_adj_ib]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                pos_adj_ib  = block_to_proc_mapping[adj_ib][2]
                if owner_adj_ib == i_proc
                    field_name[my_ib].val[1:ngx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
                        derefine_boundary_maxx(grid, field_name[pos_adj_ib].val)
                end
            end
        end

        # xmax boundary
        bc = block.boundaries[2]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                    field_name[pos_adj_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            end
        elseif bc.ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            min_y = block.min_coord_max_ref_y
            min_z = block.min_coord_max_ref_z
            min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
            min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
            if owner_adj_ib == i_proc
                if min_y == min_y_adj
                    if min_z == min_z_adj
                        field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_minx(grid, field_name[pos_adj_ib].val, 0, 0)
                    else
                        field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_minx(grid, field_name[pos_adj_ib].val, 0, 1)
                    end
                else
                    if min_z == min_z_adj
                        field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_minx(grid, field_name[pos_adj_ib].val, 1, 0)
                    else
                        field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                            refine_boundary_minx(grid, field_name[pos_adj_ib].val, 1, 1)
                    end
                end
            end
        else # ref_change == -1
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                iy = [0, 1, 0, 1][i_adj_ib]
                iz = [0, 0, 1, 1][i_adj_ib]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                pos_adj_ib  = block_to_proc_mapping[adj_ib][2]
                if owner_adj_ib == i_proc
                    field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
                        derefine_boundary_minx(grid, field_name[pos_adj_ib].val)
                end
            end
        end
    end
    println("GHOST7 done")

    #-----------------------------------------------------------------------------
    # [8] Wait buffer filling from other procs
    #-----------------------------------------------------------------------------
    stats = MPI.Waitall!(requests)
    println("GHOST8 done")

    #-----------------------------------------------------------------------------
    # [9] Copying back from buffers to fields
    #-----------------------------------------------------------------------------
    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_ref0)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_maxx_ref0:$blocks_per_proc_recv_maxx_ref0")
        for (pos_block, block) in enumerate(blocks_per_proc_recv_maxx_ref0[pos_proc])
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[nx+ngx+1:nx+2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                buffers_recv_maxx_ref0[pos_proc][:, :, :, pos_block]
        end
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refp1)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_maxx_refp1:$blocks_per_proc_recv_maxx_refp1")
        for (pos_block, block) in enumerate(blocks_per_proc_recv_maxx_refp1[pos_proc])
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[nx+ngx+1:nx+2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                buffers_recv_maxx_refp1[pos_proc][:, :, :, pos_block]
        end
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_refm1)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_maxx_refm1:$blocks_per_proc_recv_maxx_refm1")
        for (pos_block, block_tuple) in enumerate(blocks_per_proc_recv_maxx_refm1[pos_proc])
            block = block_tuple[1]
            pos_1_4 = block_tuple[2]
            iy = [0, 1, 0, 1][pos_1_4]
            iz = [0, 0, 1, 1][pos_1_4]
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
                buffers_recv_maxx_refm1[pos_proc][:, :, :, pos_block]
        end
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_minx_ref0:$blocks_per_proc_recv_minx_ref0")
        for (pos_block, block) in enumerate(blocks_per_proc_recv_minx_ref0[pos_proc])
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                buffers_recv_minx_ref0[pos_proc][:, :, :, pos_block]
        end
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refp1)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_minx_refp1:$blocks_per_proc_recv_minx_refp1")
        for (pos_block, block) in enumerate(blocks_per_proc_recv_minx_refp1[pos_proc])
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                buffers_recv_minx_refp1[pos_proc][:, :, :, pos_block]
        end
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_refm1)
        println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, blocks_per_proc_recv_minx_refm1:$blocks_per_proc_recv_minx_refm1")
        for (pos_block, block_tuple) in enumerate(blocks_per_proc_recv_minx_refm1[pos_proc])
            block = block_tuple[1]
            pos_1_4 = block_tuple[2]
            iy = [0, 1, 0, 1][pos_1_4]
            iz = [0, 0, 1, 1][pos_1_4]
            println("KAI filling i_proc=$i_proc, pos_proc=$pos_proc, proc=$proc, block=$block, pos_block=$pos_block")
            my_ib = block_to_proc_mapping[block][2]
            field_name[my_ib].val[1:ngx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
                buffers_recv_minx_refm1[pos_proc][:, :, :, pos_block]
        end
    end
    println("GHOST9 done")

    return

    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################
    ########################################################

    MPI.Finalize(); quit()

    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        tag = 100 # minx
        buffers_send_minx_ref0[:, :, :, ind_temp , pos_adj_owner_proc] = 
            field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
        request = MPI.Isend(buffers_send_minx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    for (pos_proc, proc) in enumerate(procs_adj_blocks_maxx_ref0)
        tag = 100 # minx
        request = MPI.Irecv(buffers_recv_maxx_ref0[pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end


    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        bc = block.boundaries[1]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                    field_name[pos_adj_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            else
                println("GUPX xmin boundary: block $ib send/recv to/from $adj_ib - pos-ind $pos_adj_owner_proc -- $ind_temp")
                buffers_send_minx_ref0[:, :, :, ind_temp , pos_adj_owner_proc] = 
                    field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            end
        end

        # xmax boundary
        bc = block.boundaries[2]
        if bc.ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                    field_name[pos_adj_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            else
                println("GUPX xmin boundary: block $ib send/recv to/from $adj_ib - pos-ind $pos_adj_owner_proc -- $ind_temp")
                buffers_send_minx_ref0[:, :, :, ind_temp , pos_adj_owner_proc] = 
                    field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            end
        end
    end

    stats = MPI.Waitall!(requests)

    MPI.Finalize(); quit()


    #if n_procs_adj_blocks_minx_ref0 > 0
    #    for i = 1:n_procs_adj_blocks_minx_ref0
    #        println("KET i_proc, blocks_per_proc_recv_minx_ref0: $i_proc - [$i] => $(blocks_per_proc_recv_minx_ref0[:,i:i])")
    #    end
    #end
    if n_procs_adj_blocks_maxx_ref0 > 0
        for i = 1:n_procs_adj_blocks_maxx_ref0
            println("KET i_proc, blocks_per_proc_recv_maxx_ref0: $i_proc - [$i] => $(blocks_per_proc_recv_maxx_ref0[:,i:i])")
        end
    end


    for (pos_proc, proc) in enumerate(procs_adj_blocks_minx_ref0)
        tag = 100 # minx
        n_blocks_per_proc_to_send = ind_send_minx_ref0[pos_proc]
        request = MPI.Isend(buffers_send_minx_ref0[:, :, :, 1:n_blocks_per_proc_to_send, pos_proc], proc-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    MPI.Finalize(); quit()

    # xmin boundary
    println("xmin ghost update")
    bc = block.boundaries[1]
    println("bc, bc.ref_change: $bc, $(bc.ref_change)")
    ref_change = bc.ref_change
    if ref_change == 0
        adj_ib = bc.adj_ibs[1]
        owner_adj_ib = block_to_proc_mapping[adj_ib][1]
        pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
        if owner_adj_ib == i_proc
            println("FILLING HERE")
            field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            field_name[pos_adj_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz]
        else
            pos_adj_owner_proc = findfirst(procs_adj_blocks_minx_ref0, owner_adj_ib)
            ind_send_minx_ref0[pos_adj_owner_proc] += 1 
            ind_temp = ind_send_minx_ref0[pos_adj_owner_proc]
            println("GUPX xmin boundary: block $ib send/recv to/from $adj_ib. pos-ind $pos_adj_owner_proc -- $ind_temp")
            buffers_send_minx_ref0[:, :, :, ind_temp , pos_adj_owner_proc] = 
            field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
        end
    end

    MPI.Finalize(); quit()

    max_n_adj_procs = n_procs # how many are really they?
    safety_factor = 2.0
    max_blocks_per_proc = safety_factor * n_blocks / n_procs
    buffers_send_xmin = zeros(max_n_adj_procs, max_blocks_per_proc, ngx, ny, nz)
    buffers_send_xmin_proc_mapping = zeros(max_n_adj_procs)
    buffers_send_xmin_proc_block_mapping = zeros(max_n_adj_procs, max_blocks_per_proc)

    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        println("xmin ghost update")
        bc = block.boundaries[1]
        println("bc, bc.ref_change: $bc, $(bc.ref_change)")
        ref_change = bc.ref_change
        if ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                println("FILLING HERE")
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                field_name[pos_adj_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            else
                println("GUPX xmin boundary: block $ib send/recv to/from $adj_ib. WHY HERE")
                tag = ib # tag is the index of block which sends the values
                #     if !in(owner_adj_ib, buffers_send_xmin_proc_mapping)
                #         push!(
                #     end
                #     i_send_to_owner_adj_ib = buffers_send_xmin_indexes[owner_adj_ib]
                #     buffers_send_xmin[owner_adj_ib, i_send_to_owner_adj_ib, :, :, :] = 
                #         field_name[my_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            end
        elseif ref_change == 1
            #    adj_ib = bc.adj_ibs[1]
            #    owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            #    pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            #    min_y = block.min_coord_max_ref_y
            #    min_z = block.min_coord_max_ref_z
            #    min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
            #    min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
            #    if owner_adj_ib == i_proc
            #        println("FILLING HERE")
            #        if min_y == min_y_adj
            #            if min_z == min_z_adj
            #                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 0, 0)
            #            else
            #                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 0, 1)
            #            end
            #        else
            #            if min_z == min_z_adj
            #                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 1, 0)
            #            else
            #                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_maxx(grid, field_name[pos_adj_ib].val, 1, 1)
            #            end
            #        end
            #    else
            #        tag = ib # tag is the index of block which sends the values
            #        derefined_temp = derefine_boundary_minx(grid, field_name[my_ib].val)
            #        request = MPI.Isend(derefined_temp, owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        push!(requests, request)

            #        tag = adj_ib # tag is the index of block which sends the values
            #        buffer_xmin = zeros(ngx, ny, nz)
            #        push!(buffers_xmin, buffer_xmin)
            #        request = MPI.Irecv!(buffers_xmin[end], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        #request = MPI.Irecv!(field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        push!(requests, request)
            #    end
            #else # ref_change == -1
            #    for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
            #        iy = [0, 1, 0, 1][i_adj_ib]
            #        iz = [0, 0, 1, 1][i_adj_ib]
            #        owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            #        pos_adj_ib  = block_to_proc_mapping[adj_ib][2]
            #        if owner_adj_ib == i_proc
            #            println("FILLING HERE")
            #            field_name[my_ib].val[1:ngx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
            #                derefine_boundary_maxx(grid, field_name[pos_adj_ib].val)
            #        else
            #            println("WHY HERE")
            #            tag = ib # tag is the index of block which sends the values
            #            refined_temp = refine_boundary_minx(grid, field_name[my_ib].val, iy, iz)
            #            request = MPI.Isend(refined_temp, owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            push!(requests, request)

            #            tag = adj_ib # tag is the index of block which sends the values
            #            buffer_xmin = zeros(ngx, div(ny,2), div(nz,2))
            #            push!(buffers_xmin, buffer_xmin)
            #            request = MPI.Irecv!(buffers_xmin[end], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            #request = MPI.Irecv!(field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            push!(requests, request)
            #        end
            #    end
        end

        # xmax boundary
        bc = block.boundaries[2]
        println("bc, bc.ref_change: $bc, $(bc.ref_change)")
        ref_change = bc.ref_change
        if ref_change == 0
            println("xmax ghost update")
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            if owner_adj_ib == i_proc
                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
                field_name[pos_adj_ib].val[ngx+1:2*ngx,ngy+1:ngy+ny,ngz+1:ngz+nz]
            else
                println("GUPX xmax boundary: block $ib send/recv to/from $adj_ib")
                tag = ib # tag is the index of block which sends the values
                request = MPI.Isend(field_name[my_ib].val[nx+1:ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
                push!(requests, request)

                tag = adj_ib # tag is the index of block which sends the values
                buffer_xmax = zeros(ngx, ny, nz)
                push!(buffers_xmax, buffer_xmax)
                request = MPI.Irecv!(buffers_xmax[end], owner_adj_ib-1, tag, MPI.COMM_WORLD)
                #request = MPI.Irecv!(field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
                push!(requests, request)
            end
        elseif ref_change == 1
            #    adj_ib = bc.adj_ibs[1]
            #    owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            #    pos_adj_ib   = block_to_proc_mapping[adj_ib][2]
            #    min_y = block.min_coord_max_ref_y
            #    min_z = block.min_coord_max_ref_z
            #    min_y_adj = grid.blocks[adj_ib].min_coord_max_ref_y
            #    min_z_adj = grid.blocks[adj_ib].min_coord_max_ref_z
            #    if owner_adj_ib == i_proc
            #        println("FILLING HERE")
            #        if min_y == min_y_adj
            #            if min_z == min_z_adj
            #                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_minx(grid, field_name[pos_adj_ib].val, 0, 0)
            #            else
            #                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_minx(grid, field_name[pos_adj_ib].val, 0, 1)
            #            end
            #        else
            #            if min_z == min_z_adj
            #                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_minx(grid, field_name[pos_adj_ib].val, 1, 0)
            #            else
            #                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = 
            #                    refine_boundary_minx(grid, field_name[pos_adj_ib].val, 1, 1)
            #            end
            #        end
            #    else
            #        tag = ib # tag is the index of block which sends the values
            #        derefined_temp = derefine_boundary_maxx(grid, field_name[my_ib].val)
            #        request = MPI.Isend(derefined_temp, owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        push!(requests, request)

            #        tag = adj_ib # tag is the index of block which sends the values
            #        buffer_xmax = zeros(ngx, ny, nz)
            #        push!(buffers_xmax, buffer_xmax)
            #        request = MPI.Irecv!(buffers_xmax[end], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        #request = MPI.Irecv!(field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #        push!(requests, request)
            #    end
            #else # ref_change == -1
            #    for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
            #        iy = [0, 1, 0, 1][i_adj_ib]
            #        iz = [0, 0, 1, 1][i_adj_ib]
            #        owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            #        pos_adj_ib  = block_to_proc_mapping[adj_ib][2]
            #        if owner_adj_ib == i_proc
            #            println("FILLING HERE")
            #            field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = 
            #                derefine_boundary_minx(grid, field_name[pos_adj_ib].val)
            #        else
            #            println("WHY HERE")
            #            tag = ib # tag is the index of block which sends the values
            #            refined_temp = refine_boundary_maxx(grid, field_name[my_ib].val, iy, iz)
            #            request = MPI.Isend(refined_temp, owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            push!(requests, request)

            #            tag = adj_ib # tag is the index of block which sends the values
            #            buffer_xmax = zeros(ngx, div(ny,2), div(nz,2))
            #            push!(buffers_xmax, buffer_xmax)
            #            request = MPI.Irecv!(buffers_xmax[end], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            #request = MPI.Irecv!(field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz], owner_adj_ib-1, tag, MPI.COMM_WORLD)
            #            push!(requests, request)
            #        end
            #    end
        end
    end

    for (adj_ib, buffer_to_send) in buffers_send_xmin
        block = grid.blocks[ib]
        owner_adj_ib = block_to_proc_mapping[adj_ib][1]
        request = MPI.Isend(buffer_to_send, owner_adj_ib-1, tag, MPI.COMM_WORLD)
        push!(requests, request)
    end

    stats = MPI.Waitall!(requests)

    for (ib, my_ib) in my_blocks
        block = grid.blocks[ib]

        # xmin boundary
        bc = block.boundaries[1]
        println("bc, bc.ref_change: $bc, $(bc.ref_change)")
        ref_change = bc.ref_change
        if ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib == i_proc
            else
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = splice!(buffers_xmin,1)
            end
        elseif ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib == i_proc
            else
                field_name[my_ib].val[1:ngx,ngy+1:ngy+ny,ngz+1:ngz+nz] = splice!(buffers_xmin,1)
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                iy = [0, 1, 0, 1][i_adj_ib]
                iz = [0, 0, 1, 1][i_adj_ib]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                pos_adj_ib  = block_to_proc_mapping[adj_ib][2] # useless alwas pos_ maybe
                if owner_adj_ib == i_proc
                else
                    field_name[my_ib].val[1:ngx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = splice!(buffers_xmin,1)
                end
            end
        end

        # xmax boundary
        bc = block.boundaries[2]
        println("bc, bc.ref_change: $bc, $(bc.ref_change)")
        ref_change = bc.ref_change
        if ref_change == 0
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib == i_proc
            else
                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = splice!(buffers_xmax,1)
            end
        elseif ref_change == 1
            adj_ib = bc.adj_ibs[1]
            owner_adj_ib = block_to_proc_mapping[adj_ib][1]
            if owner_adj_ib == i_proc
            else
                field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1:ngy+ny,ngz+1:ngz+nz] = splice!(buffers_xmax,1)
            end
        else
            for (i_adj_ib, adj_ib) in enumerate(bc.adj_ibs)
                iy = [0, 1, 0, 1][i_adj_ib]
                iz = [0, 0, 1, 1][i_adj_ib]
                owner_adj_ib = block_to_proc_mapping[adj_ib][1]
                pos_adj_ib  = block_to_proc_mapping[adj_ib][2] # useless alwas pos_ maybe
                if owner_adj_ib == i_proc
                else
                    field_name[my_ib].val[ngx+nx+1:2*ngx+nx,ngy+1+iy*div(ny,2):ngy+(iy+1)*div(ny,2),ngz+1+iz*div(nz,2):ngz+(iz+1)*div(nz,2)] = splice!(buffers_xmax,1)
                end
            end
        end
    end

end
