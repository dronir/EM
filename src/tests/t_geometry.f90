program t_geometry
  use geometry

  implicit none

  type(ray)                   :: r
  type(prt_sphere)            :: s
  type(intersection_geometry) :: g

  real(fd)                    :: v(3)

  call ray_init(r, RAY_TYPE_ENERGY)

  v = [1.0, 2.0, 3.0]

  print *, vec_length(v),     sqrt(14.0)
  print *, vec_length_sqr(v), 14.0

  call vec_normalize(v)

  print *, vec_length(v),     1.0
  print *, vec_length_sqr(v), 1.0  

end program t_geometry
