   subroutine SAVE_XH5F_FIELD_3D_KKP(self, xh5f, block_name, ijk, nijk, q, q_name)
   !< Save q-vector/s-scalar fields by XH5F file handler, rank 3, kind KKP.
   class(adam_object),     intent(inout)        :: self                !< ADAM.
   type(xh5f_file_object), intent(inout)        :: xh5f                !< XH5F file handler.
   character(*),           intent(in)           :: block_name          !< Block name.
   integer(I4P),           intent(in)           :: ijk(2,3)            !< Blocks extents.
   integer(I8P),           intent(in)           :: nijk(3)             !< Blocks dimensions.
   VARTYPE(KKP),           intent(in)           :: q(1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:) !< Scalar field [ni,nj,nk].
   character(*),           intent(in), optional :: q_name              !< Scalar field name.
   character(:),           allocatable          :: q_name_             !< Scalar field name, local var.
          
   q_name_ = 'scalar' ; if (present(q_name)) q_name_ = trim(adjustl(q_name))
   call xh5f%save_block_field(xdmf_field_name = q_name_,                                                    &
                              nd              = nijk,                                                       &
                              field           = q(1:ijk(1,1):ijk(2,1),ijk(1,2):ijk(2,2),ijk(1,3):ijk(2,3)), &
                              field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,                      &
                              field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,            &
                              hdf5_field_name = block_name//'-'//q_name_)
   endsubroutine SAVE_XH5F_FIELD_3D_KKP 

   subroutine SAVE_XH5F_FIELD_4D_KKP(self, xh5f, block_name, ijk, nijk, q, q_name)
   !< Save q-vector/s-scalar fields by XH5F file handler, rank 4, kind KKP.
   class(adam_object),     intent(inout)        :: self                !< ADAM.
   type(xh5f_file_object), intent(inout)        :: xh5f                !< XH5F file handler.
   character(*),           intent(in)           :: block_name          !< Block name.
   integer(I4P),           intent(in)           :: ijk(2,3)            !< Blocks extents.
   integer(I8P),           intent(in)           :: nijk(3)             !< Blocks dimensions.
   VARTYPE(KKP),           intent(in)           :: q(1:,              &
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:,&
                                                     1-self%grid%ngc:) !< Vector fields [nv,ni,nj,nk].
   character(*),           intent(in), optional :: q_name(:)           !< Vector fields names [nv].
   character(:),           allocatable          :: q_name_(:)          !< Vector fields names, local var.
   integer(I4P)                                 :: v                   !< Counter.
          
   if (present(q_name)) then
      allocate(character(len(q_name(1))):: q_name_(size(q, dim=1)))
      q_name_ = q_name
   else
      allocate(character(4):: q_name_(size(q, dim=1)))
      do v=1, size(q, dim=1)
         q_name_(v) = 'q-'//trim(strz(v,2))
      enddo
   endif
   do v=1, size(q, dim=1)
      call xh5f%save_block_field(xdmf_field_name = trim(adjustl(q_name_(v))),                                  &
                                 nd              = nijk,                                                       &
                                 field           = q(v,ijk(1,1):ijk(2,1),ijk(1,2):ijk(2,2),ijk(1,3):ijk(2,3)), &
                                 field_center    = XDMF_PARAMETERS%XDMF_ATTR_CENTER_CELL,                      &
                                 field_format    = XDMF_PARAMETERS%XDMF_DATAITEM_NUMBER_FORMAT_HDF,            &
                                 hdf5_field_name = block_name//'-'//trim(adjustl(q_name_(v))))
   enddo
   endsubroutine SAVE_XH5F_FIELD_4D_KKP 
