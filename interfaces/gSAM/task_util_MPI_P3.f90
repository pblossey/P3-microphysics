
        subroutine task_bcast_fivedim_array_real8(rank_from,array,n1,n2,n3,n4,n5)
        use mpi_stuff, only : comm
        implicit none
        include 'mpif.h'

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4,n5) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4,n5        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i1,i2,i3,i4,i5
        real :: tmp1, tmp2

        nn = n1*n2*n3*n4*n5
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_fivedim_array_real8: alloc rtmp failed!")

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do i5 = 1,n5
            do i4 = 1,n4
              do i3 = 1,n3
                do i2 = 1,n2
                  do i1 = 1,n1
                    rtmp(count) = array(i1,i2,i3,i4,i5)
                    count = count + 1
                  end do
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do i5 = 1,n5
            do i4 = 1,n4
              do i3 = 1,n3
                do i2 = 1,n2
                  do i1 = 1,n1
                    array(i1,i2,i3,i4,i5) = rtmp(count)
                    count = count + 1
                  end do
                end do
              end do
            end do
          end do
        end if

        deallocate(rtmp,STAT=ierr)
        if(ierr.gt.0) call task_abort_msg("task_bcast_fivedim_array_real8: dealloc rtmp failed!")

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,993) myrank, tmp1, tmp2
!!$        993 format('Consistency check in task_bcast_fivedim: rank/sum/sum2 = ',I4,2E16.8)

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_sixdim_array_real8(rank_from,array,n1,n2,n3,n4,n5,n6)
        use mpi_stuff, only : comm
        implicit none
        include 'mpif.h'

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4,n5,n6) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4,n5,n6        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i1,i2,i3,i4,i5,i6
        real(8) :: tmp1, tmp2

        nn = n1*n2*n3*n4*n5*n6
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_sixdim_array_real8: alloc rtmp failed!")

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
          count = 1
          do i6 = 1,n6
            do i5 = 1,n5
              do i4 = 1,n4
                do i3 = 1,n3
                  do i2 = 1,n2
                    do i1 = 1,n1
                      rtmp(count) = array(i1,i2,i3,i4,i5,i6)
                      count = count + 1
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
          count = 1
          do i6 = 1,n6
            do i5 = 1,n5
              do i4 = 1,n4
                do i3 = 1,n3
                  do i2 = 1,n2
                    do i1 = 1,n1
                      array(i1,i2,i3,i4,i5,i6) = rtmp(count)
                      count = count + 1
                    end do
                  end do
                end do
              end do
            end do
          end do
        end if

        deallocate(rtmp,STAT=ierr)
        if(ierr.gt.0) call task_abort_msg("task_bcast_sixdim_array_real8: dealloc rtmp failed!")

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,994) myrank, tmp1, tmp2
!!$        994 format('Consistency check in task_bcast_sixdim: rank/sum/sum2 = ',I4,2E16.8)

        return
        end

!----------------------------------------------------------------------

        subroutine task_bcast_sevendim_array_real8(rank_from,array,n1,n2,n3,n4,n5,n6,n7)
        use mpi_stuff, only : comm
        implicit none
        include 'mpif.h'

        integer, intent(in) :: rank_from          ! broadcasting task's rank
        real(8), intent(inout) :: array(n1,n2,n3,n4,n5,n6,n7) ! array to be broadcast
        integer, intent(in) :: n1,n2,n3,n4,n5,n6,n7        ! dimension lengths
        integer ierr

        real(8), allocatable :: rtmp(:)
        integer :: count, nn, myrank, i1,i2,i3,i4,i5,i6,i7
        real(8) :: tmp1, tmp2

        nn = n1*n2*n3*n4*n5*n6*n7
        allocate(rtmp(nn),STAT=ierr); if(ierr.gt.0) call task_abort_msg("task_bcast_sevendim_array_real8: alloc rtmp failed!")

        call MPI_COMM_RANK(comm,myrank,ierr)

        if(myrank.eq.rank_from) then
           count = 1
           do i7 = 1,n7
              do i6 = 1,n6
                 do i5 = 1,n5
                    do i4 = 1,n4
                       do i3 = 1,n3
                          do i2 = 1,n2
                             do i1 = 1,n1
                                rtmp(count) = array(i1,i2,i3,i4,i5,i6,i7)
                                count = count + 1
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end if

        call task_bcast_real8(rank_from,rtmp,nn)

        if(myrank.ne.rank_from) then
           count = 1
           do i7 = 1,n7
              do i6 = 1,n6
                 do i5 = 1,n5
                    do i4 = 1,n4
                       do i3 = 1,n3
                          do i2 = 1,n2
                             do i1 = 1,n1
                                array(i1,i2,i3,i4,i5,i6,i7) = rtmp(count)
                                count = count + 1
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end if

        deallocate(rtmp,STAT=ierr)
        if(ierr.gt.0) call task_abort_msg("task_bcast_sevendim_array_real8: dealloc rtmp failed!")

!!$        tmp1 = SUM(rtmp(:))
!!$        tmp2 = SUM(rtmp(:)*rtmp(:))
!!$        write(*,994) myrank, tmp1, tmp2
!!$        994 format('Consistency check in task_bcast_sevendim: rank/sum/sum2 = ',I4,2E16.8)

        return
        end


