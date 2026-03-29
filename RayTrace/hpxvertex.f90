module subhealpix
  use iso_c_binding
  implicit none
  contains
  subroutine getHealPixVertex(c_nside, c_ipnest, c_vertex) bind(C, name="getVertex")
    use healpix_types
    use fitstools
    use misc_utils
    use pix_tools
    implicit none
    integer(c_int), value :: c_nside
    integer(c_int), value :: c_ipnest
    real(c_float), intent(out) :: c_vertex(3,4)
    double precision :: vec(3)
    double precision :: vertex(3,4)
    integer ::i,j, nside, ipnest

    ipnest = c_ipnest
    nside = c_nside


    call pix2vec_nest(nside, ipnest, vec, vertex)

    do i = 1, 4
    do j = 1, 3
       c_vertex(j,i) = vertex(j,i)
    enddo
    enddo
  end subroutine getHealPixVertex
end module subhealpix
