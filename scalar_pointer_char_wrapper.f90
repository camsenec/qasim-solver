module scalar_pointer_char_wrapper_m
  Use iso_c_binding
  Implicit None
  Private
  Public c_charptr_to_f_charptr
  !
  ! Utility routine for getting character pointers from C.
  !
Contains
  !
  ! Utility routine to turn a C pointer to a null-terminated string
  ! into a Fortran CHARACTER pointer to that string.  The function
  ! returns a deferred-length CHARACTER pointer that is associated with
  ! the C string, and whose length (LEN) is the length of the string.
  !
  function c_charptr_to_f_charptr(ccp) Result(result)
    Type(C_ptr),Intent(In),Value :: ccp
    Character(:,C_char),Pointer :: result
    Interface
      Function strlen(p) Bind(C)
        Import C_ptr,C_size_t
        Type(C_ptr),Value :: p
        Integer(C_size_t) strlen
      End Function
    End Interface
    result => convert_cptr(ccp,strlen(ccp))
  Contains
    !
    ! This uses a variable-length CHARACTER pointer because the
    ! function C_F_pointer has no other way of encoding the length.
    !
    Function convert_cptr(p,len)
      Type(C_ptr),Intent(In) :: p
      Integer(C_size_t),Intent(In) :: len
      Character(len,C_char),Pointer :: convert_cptr
      Call C_F_pointer(p,convert_cptr)
    End Function
  End Function

end module scalar_pointer_char_wrapper_m
