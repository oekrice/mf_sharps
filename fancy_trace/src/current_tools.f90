!*******************************************************************************
MODULE current_tools
!*******************************************************************************
! Module for calculating and interpolating the current density, properly using the staggered grid
!*******************************************************************************
    USE shared_data
    USE netcdf

    IMPLICIT NONE

    contains

    subroutine calculate_current()
        !Using the read-in magnetic field, calculates the current globally.
        !Only should be done once so speed isn't really an issue

        IMPLICIT NONE

        jx = 0.0_num; jy = 0.0_num; jz = 0.0_num

        jx(0:nx+1, 0:ny,0:nz) = (bz(0:nx+1,1:ny+1,0:nz) - bz(0:nx+1, 0:ny,0:nz))/dy - (by(0:nx+1,0:ny,1:nz+1) - by(0:nx+1,0:ny,0:nz))/dz

        jy(0:nx, 0:ny+1,0:nz) = (bx(0:nx,0:ny+1,1:nz+1) - bx(0:nx, 0:ny+1,0:nz))/dz - (bz(1:nx+1,0:ny+1,0:nz) - bz(0:nx,0:ny+1,0:nz))/dx

        jz(0:nx, 0:ny,0:nz+1) = (by(1:nx+1,0:ny,0:nz+1) - by(0:nx, 0:ny,0:nz+1))/dx - (bx(0:nx,1:ny+1,0:nz+1) - bx(0:nx,0:ny,0:nz+1))/dy

        print*, jx(nx/2,ny/2,0:nz+1)
        print*, maxval(abs(jx(0:nx+1,0:ny,2:nz-2))), maxval(abs(jy)), maxval(abs(jz))
        read(*,*)

    end subroutine calculate_current

    SUBROUTINE interpolate_jfield(x_point, y_point, z_point)
    !Finds the square of the current field at a given point, properly using the staggered grid as per
    IMPLICIT NONE

    REAL(num):: x_point, y_point, z_point !Coordinates of point to interpolate
    REAL(num):: xp, yp, zp !Coordinates in the respective dimentions
    INTEGER:: xi, yi, zi !Cell indices in each dimension
    REAL(num):: xf, yf, zf !Distance 'up' each cell

    j1 = 0.0_num
    !Establish ratios
    xp = nx*(x_point - x0)/(x1 - x0)
    yp = ny*(y_point - y0)/(y1 - y0)
    zp = nz*(z_point - z0)/(z1 - z0)

    !Interpolate jx
    xi = int(xp + 0.5_num); yi = int(yp); zi = int(zp)
    xf = xp + 0.5_num - xi; yf = yp - yi; zf = zp - zi

    j1(0) = j1(0) + jx(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + jx(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    j1(0) = j1(0) + jx(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + jx(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    j1(0) = j1(0) + jx(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + jx(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    j1(0) = j1(0) + jx(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + jx(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate jy
    xi = int(xp); yi = int(yp+0.5_num); zi = int(zp)
    xf = xp - xi; yf = yp + 0.5_num - yi; zf = zp  - zi

    j1(1) = j1(1) + jy(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + jy(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    j1(1) = j1(1) + jy(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + jy(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    j1(1) = j1(1) + jy(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + jy(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    j1(1) = j1(1) + jy(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + jy(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    !Interpolate jz
    xi = int(xp); yi = int(yp); zi = int(zp + 0.5_num)
    xf = xp - xi; yf = yp - yi; zf = zp + 0.5_num - zi

    j1(2) = j1(2) + jz(xi,yi,zi)*(1.0_num-xf)*(1.0_num-yf)*(1.0_num-zf) + jz(xi,yi,zi+1)*(1.0_num-xf)*(1.0_num-yf)*(zf)
    j1(2) = j1(2) + jz(xi,yi+1,zi)*(1.0_num-xf)*(yf)*(1.0_num-zf)       + jz(xi,yi+1,zi+1)*(1.0_num-xf)*(yf)*(zf)
    j1(2) = j1(2) + jz(xi+1,yi,zi)*(xf)*(1.0_num-yf)*(1.0_num-zf)       + jz(xi+1,yi,zi+1)*(xf)*(1.0_num-yf)*(zf)
    j1(2) = j1(2) + jz(xi+1,yi+1,zi)*(xf)*(yf)*(1.0_num-zf)             + jz(xi+1,yi+1,zi+1)*(xf)*(yf)*(zf)

    jmag = sum(j1**2)

    END SUBROUTINE interpolate_jfield

END MODULE current_tools
