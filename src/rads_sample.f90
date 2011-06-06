use rads                    ! Use rads module
type(rads_sat) :: S         ! Define a structure that contains all mission information
type(rads_pass) :: P        ! Structure for pass information
real(eightbytereal), allocatable :: lon(:), lat(:), sla(:)
integer :: cycle, pass, i
call rads_init (S, 'e1')    ! Initialize RADS. Select ERS-1. No mission needed!
do cycle = 1,30
    do pass = 1,10
        call rads_open_pass (S, P, cycle, pass)              ! Open pass file
        allocate (lon(P%ndata), lat(P%ndata), sla(P%ndata))  ! Allocate memory
        call rads_get_var (S, P, 'lon', lon)                 ! Get the data
        call rads_get_var (S, P, 'lat', lat)
        call rads_get_var (S, P, 'sla', sla)
        do i = 1,P%ndata
            write (*,*) lon(i), lat(i), sla(i)               ! Print the data
        enddo
        deallocate (lon, lat, sla)                           ! Deallocate memory
        call rads_close_pass (S, P)                          ! Close pass
    enddo
enddo
call rads_end (S)                                            ! End
end
