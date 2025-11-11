   subroutine REGISTER_AUX_FIELD_4D_KKP(q_src, q_name_src, q_reg, q_name_reg)
   !< Register auxiliary fields data into ADAM IO class, rank 4, kind KKP.
   VARTYPE(KKP), intent(in), target   :: q_src(:,:,:,:) !< Auxiliary vector cell centered fields, source.
   character(*), intent(in), optional :: q_name_src     !< Auxiliary vector fields names, source.
   VARTYPE(KKP), intent(out), pointer :: q_reg(:,:,:,:) !< Auxiliary vector cell centered fields, registered.
   type(string), intent(out)          :: q_name_reg     !< Auxiliary vector cell centered fields names, registered.

   q_reg => null()
   q_reg => q_src
   if (present(q_name_src)) then
      q_name_reg = trim(adjustl(q_name_src))
   else
      q_name_reg = 'scalar'
   endif
   endsubroutine REGISTER_AUX_FIELD_4D_KKP

   subroutine REGISTER_AUX_FIELD_5D_KKP(q_src, q_name_src, q_reg, q_name_reg)
   !< Register auxiliary fields data into ADAM IO class, rank 5, kind KKP.
   VARTYPE(KKP), intent(in), target       :: q_src(:,:,:,:,:) !< Auxiliary vector cell centered fields, source.
   character(*), intent(in), optional     :: q_name_src(:)    !< Auxiliary vector fields names, source.
   VARTYPE(KKP), intent(out), pointer     :: q_reg(:,:,:,:,:) !< Auxiliary vector cell centered fields, registered.
   type(string), intent(out), allocatable :: q_name_reg(:)    !< Auxiliary vector cell centered fields names, registered.
   integer(I4P)                           :: v                !< Counter.

   q_reg => null()
   q_reg => q_src
   allocate(q_name_reg(size(q_src, dim=1)))
   if (present(q_name_src)) then
      do v=1, size(q_src, dim=1)
         q_name_reg(v) = trim(adjustl(q_name_src(v)))
      enddo
   else
      do v=1, size(q_src, dim=1)
         q_name_reg(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   endsubroutine REGISTER_AUX_FIELD_5D_KKP

   subroutine SAVE_XH5F_FIELD_4D_KKP(self, xh5f, block_name, ijk, nijk, q, q_name)
   !< Save q-vector/s-scalar fields by XH5F file handler, rank 4, kind KKP.
   class(io_object),       intent(inout) :: self                !< IO handler.
   type(xh5f_file_object), intent(inout) :: xh5f                !< XH5F file handler.
   character(*),           intent(in)    :: block_name          !< Block name.
   integer(I4P),           intent(in)    :: ijk(2,3)            !< Blocks extents.
   integer(I8P),           intent(in)    :: nijk(3)             !< Blocks dimensions.
   VARTYPE(KKP),           intent(in)    :: q(1-self%grid%ngc:,&
                                              1-self%grid%ngc:,&
                                              1-self%grid%ngc:) !< Scalar field [ni,nj,nk].
   type(string),           intent(in)    :: q_name              !< Scalar field name.
          
   call xh5f%save_block_field(xdmf_field_name = q_name%chars(),                                           &
                              nd              = nijk,                                                     &
                              field           = q(ijk(1,1):ijk(2,1),ijk(1,2):ijk(2,2),ijk(1,3):ijk(2,3)), &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,                    &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,          &
                              hdf5_field_name = block_name//'-'//q_name%chars())
   endsubroutine SAVE_XH5F_FIELD_4D_KKP 

   subroutine SAVE_XH5F_FIELD_5D_KKP(self, xh5f, block_name, ijk, nijk, q, q_name)
   !< Save q-vector/s-scalar fields by XH5F file handler, rank 5, kind KKP.
   class(io_object),       intent(inout) :: self                !< IO handler.
   type(xh5f_file_object), intent(inout) :: xh5f                !< XH5F file handler.
   character(*),           intent(in)    :: block_name          !< Block name.
   integer(I4P),           intent(in)    :: ijk(2,3)            !< Blocks extents.
   integer(I8P),           intent(in)    :: nijk(3)             !< Blocks dimensions.
   VARTYPE(KKP),           intent(in)    :: q(1:,              &
                                              1-self%grid%ngc:,&
                                              1-self%grid%ngc:,&
                                              1-self%grid%ngc:) !< Vector fields [nv,ni,nj,nk].
   type(string),           intent(in)    :: q_name(:)           !< Vector fields names [nv].
   integer(I4P)                          :: v                   !< Counter.
          
   do v=1, size(q, dim=1)
      call xh5f%save_block_field(xdmf_field_name = q_name(v)%chars(),                                          &
                                 nd              = nijk,                                                       &
                                 field           = q(v,ijk(1,1):ijk(2,1),ijk(1,2):ijk(2,2),ijk(1,3):ijk(2,3)), &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,                      &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,            &
                                 hdf5_field_name = block_name//'-'//q_name(v)%chars())
   enddo
   endsubroutine SAVE_XH5F_FIELD_5D_KKP 
