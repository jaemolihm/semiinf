module iter_slab
    use comms
    use parameters
    use hamiltonian, ONLY : h00, h01, h11, h12, omega
    use postprocess_green, ONLY : green_s, green_s1, green_b
    IMPLICIT NONE
    SAVE
    PRIVATE
    ! private temporary variables used in iter_slab_update and iter_slab_green
    COMPLEX(dp), allocatable, dimension(:,:) :: a_s, a_b, b_s, b_b, e_s0, e_s1, e_b
    complex(DP), allocatable, dimension(:,:) :: inv_e_b, a_b_prev, b_b_prev, &
        a_s_prev, b_s_prev, temp_mat, temp_mat_sb
    integer, allocatable, dimension(:) :: ipiv_s, ipiv_b
    ! private variables used in iter_slab_check_convergence
    REAL(dp) :: conver_as, conver_ab, conver_bs, conver_bb
    PUBLIC :: iter_slab_main, iter_slab_allocate, iter_slab_deallocate, iter_slab_green_finite
contains
    SUBROUTINE iter_slab_main()
        CALL iter_slab_initialize()
        ! ... iterate semi-infinite slab
        DO n_iter = 1, MAX_N_ITER
            CALL iter_slab_update()
            CALL iter_slab_check_convergence()
            IF (flag_converged) EXIT
        ENDDO
        CALL iter_slab_green()
    END SUBROUTINE

    SUBROUTINE iter_slab_initialize()
        ! ... initialize iteration paramters
        IMPLICIT NONE
        e_s0 = h00
        e_s1 = h11
        e_b = h11
        a_s = h01
        b_s = TRANSPOSE(CONJG(h01))
        a_b = h12
        b_b = TRANSPOSE(CONJG(h12))
    END SUBROUTINE

    SUBROUTINE iter_slab_check_convergence()
        ! check convergences of hopping matrices
        IMPLICIT NONE
        INTEGER :: i, j
        conver_ab = 0.0_dp
        conver_bb = 0.0_dp
        conver_as = 0.0_dp
        conver_bs = 0.0_dp
        flag_converged = .false.
        DO j = 1, nbulk
            DO i = 1, nbulk
                conver_ab = conver_ab + SQRT(REAL(a_b(i,j),dp)**2+AIMAG(a_b(i,j))**2)
                conver_bb = conver_bb + SQRT(REAL(b_b(i,j),dp)**2+AIMAG(b_b(i,j))**2)
            END DO
        END DO
        DO j = 1, nbulk
            DO i = 1, nsurf
                conver_as = conver_as + SQRT(REAL(a_s(i,j),dp)**2+AIMAG(a_s(i,j))**2)
            END DO
        END DO
        DO j = 1, nsurf
            DO i = 1, nbulk
                conver_bs = conver_bs + SQRT(REAL(b_s(i,j),dp)**2+AIMAG(b_s(i,j))**2)
            END DO
        END DO
        IF (conver_ab.LT.hopping_tol .AND. conver_bb.LT.hopping_tol &
            .AND. conver_as.LT.hopping_tol .AND. conver_bs.LT.hopping_tol) THEN
            flag_converged = .true.
        END IF
    END SUBROUTINE

    SUBROUTINE iter_slab_update()
        ! use recurrence relation to update matrices
        IMPLICIT NONE
        ! inv_e_b = (omgea - e_b)^-1
        CALL inv_omega_minus_mat(nbulk, e_b, omega, inv_e_b, 'iter_bulk_update for inv_e_b')

        CALL ZCOPY(nsurf*nbulk, a_s, 1, a_s_prev, 1)
        CALL ZCOPY(nbulk*nsurf, b_s, 1, b_s_prev, 1)
        CALL ZCOPY(nbulk*nbulk, a_b, 1, a_b_prev, 1)
        CALL ZCOPY(nbulk*nbulk, b_b, 1, b_b_prev, 1)

        ! temp_mat = inv_e_b * a_b
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,inv_e_b,nbulk,a_b,nbulk,czero,temp_mat,nbulk)
        ! a_s = a_s_prev * temp_mat
        CALL ZGEMM('N','N',nsurf,nbulk,nbulk,cone,a_s_prev,nsurf,temp_mat,nbulk,czero,a_s,nsurf)
        ! a_b = a_b_prev * temp_mat
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,czero,a_b,nbulk)
        ! e_b = e_b(prev) + b_b_prev * temp_mat
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,b_b_prev,nbulk,temp_mat,nbulk,cone,e_b,nbulk)

        ! temp_mat = b_b_prev * inv_e_b
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,b_b_prev,nbulk,inv_e_b,nbulk,czero,temp_mat,nbulk)
        ! b_s = temp_mat * b_s_prev
        CALL ZGEMM('N','N',nbulk,nsurf,nbulk,cone,temp_mat,nbulk,b_s_prev,nbulk,czero,b_s,nbulk)
        ! b_b = temp_mat * b_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,temp_mat,nbulk,b_b_prev,nbulk,czero,b_b,nbulk)

        ! temp_mat = inv_e_b * b_b_prev
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,inv_e_b,nbulk,b_b_prev,nbulk,czero,temp_mat,nbulk)
        ! e_b = e_b(prev) + a_b_prev * temp_mat
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,cone,e_b,nbulk)
        ! e_s1 = e_s1(prev) + a_b_prev * temp_mat
        CALL ZGEMM('N','N',nbulk,nbulk,nbulk,cone,a_b_prev,nbulk,temp_mat,nbulk,cone,e_s1,nbulk)

        ! temp_mat_sb = a_s_prev * inv_e_b
        CALL ZGEMM('N','N',nsurf,nbulk,nbulk,cone,a_s_prev,nsurf,inv_e_b,nbulk,czero,temp_mat_sb,nsurf)
        ! e_s0 = e_s0(prev) + temp_mat_sb * b_s_prev
        CALL ZGEMM('N','N',nsurf,nsurf,nbulk,cone,temp_mat_sb,nsurf,b_s_prev,nbulk,cone,e_s0,nsurf)
    END SUBROUTINE

    SUBROUTINE iter_slab_green()
        ! calculate green function after iteration
        IMPLICIT NONE
        CALL inv_omega_minus_mat(nsurf, e_s0, omega, green_s, 'iter_slab_green for green_s')
        CALL inv_omega_minus_mat(nbulk, e_s1, omega, green_s1, 'iter_slab_green for green_s1')
        CALL inv_omega_minus_mat(nbulk, e_b, omega, green_b, 'iter_slab_green for green_b')
    END SUBROUTINE

    SUBROUTINE iter_slab_green_finite()
        ! calculate green function for finite slab by direct inversion
        IMPLICIT NONE
        CALL inv_omega_minus_mat(nsurf, h00, omega, green_s, 'iter_slab_green for green_s')
    END SUBROUTINE

    SUBROUTINE iter_slab_allocate()
        IMPLICIT NONE
        INTEGER :: ierr
        ALLOCATE(a_s(nsurf, nbulk),stat=ierr)
        ALLOCATE(b_s(nbulk, nsurf),stat=ierr)
        ALLOCATE(a_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(b_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(e_s0(nsurf, nsurf),stat=ierr)
        ALLOCATE(e_s1(nbulk, nbulk),stat=ierr)
        ALLOCATE(e_b(nbulk, nbulk),stat=ierr)
        ALLOCATE(inv_e_b(nbulk,nbulk),stat=ierr)
        ALLOCATE(a_s_prev(nsurf,nbulk),stat=ierr)
        ALLOCATE(b_s_prev(nbulk,nsurf),stat=ierr)
        ALLOCATE(a_b_prev(nbulk,nbulk),stat=ierr)
        ALLOCATE(b_b_prev(nbulk,nbulk),stat=ierr)
        ALLOCATE(temp_mat(nbulk,nbulk),stat=ierr)
        ALLOCATE(temp_mat_sb(nsurf,nbulk),stat=ierr)
        ALLOCATE(ipiv_s(nsurf),stat=ierr)
        ALLOCATE(ipiv_b(nbulk),stat=ierr)
        if (ierr /= 0) CALL io_error('Error in allocation in iter_slab_allocate')
    END SUBROUTINE

    SUBROUTINE iter_slab_deallocate()
        IMPLICIT NONE
        INTEGER :: ierr
        DEALLOCATE(a_s,stat=ierr)
        DEALLOCATE(b_s,stat=ierr)
        DEALLOCATE(a_b,stat=ierr)
        DEALLOCATE(b_b,stat=ierr)
        DEALLOCATE(e_s0,stat=ierr)
        DEALLOCATE(e_s1,stat=ierr)
        DEALLOCATE(e_b,stat=ierr)
        DEALLOCATE(inv_e_b,stat=ierr)
        DEALLOCATE(a_s_prev,stat=ierr)
        DEALLOCATE(b_s_prev,stat=ierr)
        DEALLOCATE(a_b_prev,stat=ierr)
        DEALLOCATE(b_b_prev,stat=ierr)
        DEALLOCATE(temp_mat,stat=ierr)
        DEALLOCATE(temp_mat_sb,stat=ierr)
        DEALLOCATE(ipiv_s,stat=ierr)
        DEALLOCATE(ipiv_b,stat=ierr)
        if (ierr /= 0) CALL io_error('Error in deallocation in iter_slab_deallocate')
    END SUBROUTINE
end module iter_slab
