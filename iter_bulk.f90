module iter_bulk
    use constants
    use parameters
    use hamiltonian, ONLY : h11, h12, omega
    use postprocess_green, ONLY : green_s, green_b
    IMPLICIT NONE
    SAVE
    PRIVATE
    ! private temporary variables used in iter_bulk_update and iter_bulk_green
    COMPLEX(dp), allocatable, dimension(:,:) :: a_b, b_b, e_s, e_b
    complex(DP), allocatable, dimension(:,:) :: inv_e_b, a_b_prev, b_b_prev, temp_mat
    integer, allocatable, dimension(:) :: ipiv
    ! private variables used in iter_bulk_check_convergence
    REAL(dp) :: conver1, conver2
    PUBLIC :: iter_bulk_main, iter_bulk_allocate, iter_bulk_deallocate
contains
    SUBROUTINE iter_bulk_main()
        CALL iter_bulk_initialize()
        ! ... iterate semi-infinite slab
        DO n_iter = 1, MAX_N_ITER
            CALL iter_bulk_update()
            CALL iter_bulk_check_convergence()
            IF (flag_converged) EXIT
        ENDDO
        CALL iter_bulk_green()
    END SUBROUTINE

    SUBROUTINE iter_bulk_initialize()
        ! ... initialize iteration paramters
        IMPLICIT NONE
        e_s = h11
        e_b = h11
        a_b = h12
        b_b = TRANSPOSE(CONJG(h12))
    END SUBROUTINE

    SUBROUTINE iter_bulk_check_convergence()
        ! check convergences of hopping matrices
        IMPLICIT NONE
        INTEGER :: i, j
        conver1 = 0.0_dp
        conver2 = 0.0_dp
        flag_converged = .false.
        DO j = 1, nbulk
            DO i = 1, nbulk
                conver1 = conver1 + SQRT(REAL(a_b(i,j),dp)**2+AIMAG(a_b(i,j))**2)
                conver2 = conver2 + SQRT(REAL(b_b(i,j),dp)**2+AIMAG(b_b(i,j))**2)
            END DO
        END DO
        IF (conver1.LT.hopping_tol .AND. conver2.LT.hopping_tol) THEN
            flag_converged = .true.
        END IF
    END SUBROUTINE

    SUBROUTINE iter_bulk_update()
        ! use recurrence relation to update matrices
        IMPLICIT NONE
        ! inv_e_b = (omgea - e_b)^-1
        CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, 'iter_bulk_update for inv_e_b')

        CALL ZCOPY(nbulk*nbulk, a_b, 1, a_b_prev, 1)
        CALL ZCOPY(nbulk*nbulk, b_b, 1, b_b_prev, 1)

        ! temp_mat = a_b_prev * inv_e_b
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,a_b_prev,nbulk,inv_e_b,nbulk,CMPLX_0,temp_mat,nbulk)
        ! e_s = e_s(prev) + temp_mat * b_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,temp_mat,nbulk,b_b_prev,nbulk,CMPLX_1,e_s,nbulk)
        ! e_b = e_b(prev) + temp_mat * b_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,temp_mat,nbulk,b_b_prev,nbulk,CMPLX_1,e_b,nbulk)
        ! a_b = temp_mat * a_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,temp_mat,nbulk,a_b_prev,nbulk,CMPLX_0,a_b,nbulk)   

        ! temp_mat = b_b_prev * inv_e_b
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,b_b_prev,nbulk,inv_e_b,nbulk,CMPLX_0,temp_mat,nbulk)
        ! e_b = e_b(prev) + temp_mat * a_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,temp_mat,nbulk,a_b_prev,nbulk,CMPLX_1,e_b,nbulk)
        ! b_b = temp_mat * b_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,CMPLX_1,temp_mat,nbulk,b_b_prev,nbulk,CMPLX_0,b_b,nbulk)
    END SUBROUTINE

    SUBROUTINE iter_bulk_green()
        ! calculate green function after iteration
        IMPLICIT NONE
        CALL inv_omega_minus_mat(nbulk, e_s, omega, green_s, 'iter_bulk_green for green_s')
        CALL inv_omega_minus_mat(nbulk, e_b, omega, green_b, 'iter_bulk_green for green_b')
    END SUBROUTINE

    SUBROUTINE iter_bulk_allocate()
        IMPLICIT NONE
        INTEGER :: ierr
        ALLOCATE(a_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(b_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(e_s(nbulk, nbulk),stat=ierr)
        ALLOCATE(e_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(inv_e_b(nbulk,nbulk),stat=ierr)
        ALLOCATE(a_b_prev(nbulk,nbulk),stat=ierr)
        ALLOCATE(b_b_prev(nbulk,nbulk),stat=ierr)
        ALLOCATE(temp_mat(nbulk,nbulk),stat=ierr)
        ALLOCATE(ipiv(nbulk),stat=ierr)
        if (ierr /= 0) CALL io_error('Error in allocation in iter_bulk_allocate')
    END SUBROUTINE

    SUBROUTINE iter_bulk_deallocate()
        IMPLICIT NONE
        INTEGER :: ierr
        DEALLOCATE(a_b,stat=ierr)
        DEALLOCATE(b_b,stat=ierr)
        DEALLOCATE(e_s,stat=ierr)
        DEALLOCATE(e_b,stat=ierr)
        DEALLOCATE(inv_e_b,stat=ierr)
        DEALLOCATE(a_b_prev,stat=ierr)
        DEALLOCATE(b_b_prev,stat=ierr)
        DEALLOCATE(temp_mat,stat=ierr)
        DEALLOCATE(ipiv,stat=ierr)
        if (ierr /= 0) CALL io_error('Error in deallocation in iter_bulk_deallocate')
    END SUBROUTINE
end module iter_bulk
