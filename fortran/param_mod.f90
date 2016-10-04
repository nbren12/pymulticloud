module param_mod
  implicit none

  real(8) :: hour, minute, day, km, tempouttime
  parameter(hour=3600.d0, minute=60.0d0, day=86400.0d0, km=1000.d0)

  real(8) pi, cp, gamma, gammam, theta0, omega, r, g, beta, EQ, n2, zt, zm, zb, zp


contains

  subroutine init_param()
    pi=4*DATAN(1.d0)

    cp=1000.d0     !J/kg/K
    gamma=1.7D0    ! ratio of moist and dry lapse rates
    gammam = 6.0 / km ! moist lapse rate (K /m)
    theta0 = 300  ! reference pot temp (K)
    omega= 2*pi/(24*hour)! the angular velocity of the earth and
    r=6378*km           ! radius of the earth
    theta0=300.d0 !K background  reference temperature
    g=9.8D0 !m/s gravitational acceleration
    beta=2*omega/r      !1/s/m
    EQ=40000*km!  Earth's peremeter at the Equator
    n2=.0001D0         !Vaissala buoyancy frequency squared


    zt=16.00*km !topospheric height
    zm = 5000.d0 !middle trop height
    zb=500.d0 !meters boundary layer depth
    zp = 8.d0 * km ! average height of penetrative clouds
  end subroutine init_param

end module param_mod
