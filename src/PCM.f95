!    Copyright (C) 2013, Kevin Shook, Centre for Hydrology
!    University of Saskatchewan, 12 Kirk Hall, 117 Science Place, Saskatoon, SK, S7N 5C8

!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!

module convert
implicit none
real, parameter :: p1 = 1.72
real, parameter :: p2 = 3.33
real, parameter :: maxerror = 0.00001
integer, parameter :: maxiterations = 1000
contains

real function areafraction (maxarea, volume, S, h)
real :: volume, maxarea, S, h, pval, areafrac
real :: currentdepth, currentarea
! set value of scaling parameter
  if (maxarea < 10000.0) then
    pval= p1
  else
    pval = p2
  endif
  
! calculate current depth & area
  currentdepth = (volume*(1.0+2.0/pval)/S) ** (1.0/(1.0+2.0/pval))
  currentdepth = min(currentdepth, h)
  currentarea = S * (currentdepth)**(2.0/pval)
  areafrac = currentarea / maxarea
  areafrac = min(areafrac, 1.0)
  areafrac = max(areafrac, 0.0)
  areafraction = areafrac

end function areafraction

subroutine calc_h_and_S(maxvolume, maxarea, S, h)
real :: maxvolume, maxarea, S, h, est_area, area_error, pval
integer :: iterations
logical :: done
! does iterative calculations to find h and S
  if (maxarea < 10000.0) then
     pval = p1
  else
     pval = p2
  endif
  h = 1.0
  done = .false.
  iterations = 0
  do while (.not.(done))
    S = maxarea / (h ** (2.0 / pval))
    h = ((maxvolume * (1.0 + 2.0 / pval)) / S) ** (1.0 / (1.0 + 2.0/pval))
    est_area = S  * (h ** (2.0 / pval))
    area_error = abs(est_area - maxarea) / maxarea

    if ( (area_error < maxerror) .or. (iterations > maxiterations )) then
      done = .true.
    else
      iterations = iterations + 1
    endif
  end do
end subroutine calc_h_and_S

subroutine print_args
  print *, 'Copyright (C) 2012, Kevin Shook, Centre for Hydrology'
  print *, 'This program is free software: you can redistribute it and/or modify'
  print *, 'it under the terms of the GNU General Public License as published by'
  print *, 'the Free Software Foundation, either version 3 of the License, or'
  print *, '(at your option) any later version.'
  print *
  print *, 'This program is distributed in the hope that it will be useful,'
  print *, 'but WITHOUT ANY WARRANTY; without even the implied warranty of'
  print *, 'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
  print *, 'GNU General Public License for more details.'
  print *,
  print *, 'You should have received a copy of the GNU General Public License'
  print *, 'along with this program.  If not, see <http://www.gnu.org/licenses/>'

  print *, '7 arguments are required'
  print *, '1. Max area file name'
  print *, '2. Max volume file name'
  print *, '3. Connectivity file name'
  print *, '4. Initial state (full/empty/state file'
  print *, '5. Precipitation/evaporation (real)'
  print *, '6. Runoff fraction (real)'
  print *, '7. Output file name'

end subroutine print_args


real function wetland_drainage_area(maxarea)
real :: maxarea
! based on scaling relationships
 if (maxarea > 1000.0) then
  wetland_drainage_area = exp(3.44329) * (maxarea ** 0.738377)
 else
  wetland_drainage_area = 2893.887
 endif
end function wetland_drainage_area

function upcase(string) result(upper)
! copied from http://www.star.le.ac.uk/~cgp/fortran.html
character(len=*), intent(in) :: string
character(len=len(string)) :: upper
integer :: j
  do j = 1,len(string)
    if(string(j:j) >= "a" .and. string(j:j) <= "z") then
      upper(j:j) = achar(iachar(string(j:j)) - 32)
    else
      upper(j:j) = string(j:j)
    end if
  end do
end function upcase

end module convert


program event
use convert
!   This program does fill & spill calculations for a group of simulated sloughs
!   This version uses specified areas and volumes with revised scaling
!   (c) Kevin Shook, Sept. 7, 2011, Apr. 4, 2013
!   This is free-form Fortran 95
!   This program does the following:
!   1.  gets parameters from command line
!   2.  reads data from files: MaxArea, MaxVolume, Routing, Inital State
!   3.  applies precip. or evap value for a single event
!   4.  routes surplus volumes until all sloughs are not overflowing

!   April 4, 2013
!   Added willow ring area for evaporation, ROF (runoff factor) for precipitation

implicit none

character(80) :: MaxAreaFileName, RoutingFileName, OutputStateFileName, MaxVolumeFileName
character(40) :: InitialState, PrecipString, ROFString

real, allocatable, dimension(:) :: current_area, current_depth, current_volume, max_area
real, allocatable, dimension(:) :: basin_area, max_volume, max_depth, S_values, h_values, wetland_drainage
integer, allocatable, dimension(:) :: dests

real :: precip, evap_m, ROF
logical :: done
integer :: i, j, n, arg_count, destination, iter_count, NumSloughs, p, i_evap
real :: largest_value, transferred_volume, dummy
real, parameter :: tolerance=0.1  ! tolerence for convergence (m³)
real :: initial_vol, delta_vol, volfrac, areafrac, runofffrac
real :: outflow_volume, final_vol, final_area, total_basin_area, max_slough_area, max_slough_volume
integer, parameter :: iter_max = 100

arg_count = command_argument_count()
if (arg_count /= 7) then
  call print_args
  stop
end if

! extract arguments
CALL get_command_argument(1, MaxAreaFileName)
CALL get_command_argument(2, MaxVolumeFileName)
CALL get_command_argument(3, RoutingFileName)
CALL get_command_argument(4, InitialState)
CALL get_command_argument(5, PrecipString)
CALL get_command_argument(6, ROFString)
CALL get_command_argument(7, OutputStateFileName)

read(PrecipString,'(f8.1)') precip
read(ROFString,'(f8.6)') ROF

!  now open files and input arrays
!  read the draining destinations
open(1, file= RoutingFileName, access="sequential", action = "read")
n=0
do
  read(1,*,end=10) dummy
  n=n+1
enddo
10   rewind(1)

NumSloughs = n

allocate(max_area(NumSloughs))
allocate(max_volume(NumSloughs))
allocate(max_depth(NumSloughs))
allocate(dests(NumSloughs))
allocate(current_area(NumSloughs))
allocate(current_depth(NumSloughs))
allocate(current_volume(NumSloughs))
allocate(basin_area(NumSloughs))
allocate(S_values(NumSloughs))
allocate(h_values(NumSloughs))
allocate(wetland_drainage(NumSloughs))

do i =  1, NumSloughs
  read(1,*) dests(i)
end do
close(1)

!  print *, 'Reading slough max area from ', MaxAreaFileName
open(1, file= MaxAreaFileName, access="sequential", action = "read")
do i =  1, NumSloughs
  read(1,*) max_area(i)
end do
close(1)

!  print *, 'Reading slough max volume from ', MaxVolumeFileName
open(1, file= MaxVolumeFileName, access="sequential", action = "read")
do i =  1, NumSloughs
  read(1,*) max_volume(i)
end do
close(1)

max_slough_volume = sum(max_volume)
max_slough_area = sum(max_area)
do i =  1, NumSloughs
  wetland_drainage(i) = wetland_drainage_area(max_area(i))
end do

total_basin_area = sum(wetland_drainage)

! Calculate values for S and h
do i =  1, NumSloughs
  call calc_h_and_S(max_volume(i), max_area(i), S_values(i), h_values(i))
end do

! set the initial conditions
if (upcase(InitialState) == 'EMPTY') then
  current_volume = 0
elseif (upcase(InitialState) == 'FULL') then
  current_volume = max_volume
else
  ! read the draining destinations from file
  open(1, file= InitialState, access="sequential", action = "read")
  do i =  1, NumSloughs
     read(1,*) current_volume(i)
  end do
  close(1)
end if

initial_vol = sum(current_volume)

! apply precip/evap
if (precip > 0.) then
  done = .true.
  do j=1,NumSloughs
  ! apply rainfall
  ! do direct precip + runoff
  ! get current water area
    areafrac = areafraction(max_area(j), current_volume(j), S_values(j), h_values(j))

  ! do direct precip on water
    current_volume(j)=current_volume(j)+(precip*max_area(j)*areafrac)/1000.0

  ! now do runoff from unwetted upland
    current_volume(j)=current_volume(j)+&
    (precip * ROF /1000.0) * (wetland_drainage(j) - areafrac * max_area(j))

  ! print *, 'Slough ', j, 'basin area', basin_area(j)
    if (current_volume(j) >= max_volume(j)) then
      done = .false.
    end if

  end do ! adding precip, doing runoff

! route excess flows
!    print *, 'Route excess flows'
  iter_count =1
  do while (.not.(done))
    largest_value = 0.0
    do j = NumSloughs, 2, -1
      if (current_volume(j) > max_volume(j)) then
        transferred_volume = current_volume(j) - max_volume(j)
        current_volume(j) = max_volume(j)
        destination = dests(j)
        current_volume(destination) = current_volume(destination)+transferred_volume
      else
        transferred_volume = 0.0
      end if
      if (transferred_volume > largest_value) then
        largest_value=transferred_volume
      end if
    end do  ! each slough

    if ((largest_value <= tolerance) .or. (iter_count > iter_max)) then
      done = .true.
    else
      iter_count = iter_count+1
    end if
  end do ! done (all routing)
! now do slough #1
  outflow_volume = max(current_volume(1)-max_volume(1),0.0)
  current_volume(1)=min(current_volume(1), max_volume(1))
else
! do evap - evaporate water 1 mm at a time
  outflow_volume = 0.0
  i_evap = int(-1.0 * precip)
  evap_m = 1.0/1000.0
  do j=1, NumSloughs
    do p = 1, i_evap
      ! assume area of willow ring = max area of water
      current_volume(j) = max((current_volume(j) - (evap_m * max_area(j))),0.0)
    end do
  end do
end if

! calculate stats & write to stdout
final_vol = sum(current_volume)
do j=1, NumSloughs
  current_area(j) = areafraction(max_area(j), current_volume(j), S_values(j), h_values(j) ) * max_area(j)
end do

final_area = sum(current_area)
delta_vol = final_vol - initial_vol
volfrac=final_vol/max_slough_volume
areafrac=final_area/max_slough_area

! calculate fractional outflow
if (precip > 0.) then
  runofffrac = 1.0 - (delta_vol / ((precip/1000.0) * total_basin_area))
else
  runofffrac = 0.0
endif
!  print *, 'total_basin_area max_slough_area max_slough_volume', &
!            ' InitialState precip final_vol final_area mean_area delta_vol outflow_volume'
write (*,40) NumSloughs, total_basin_area, max_slough_area, max_slough_volume, InitialState, precip, &
   final_vol, final_area, delta_vol, outflow_volume, volfrac, areafrac, runofffrac
! finally, write output to file
open(3, file= OutputStateFileName, access="sequential", action = "write", status="replace")
do i = 1, NumSloughs
  write (3,20)  current_volume(i)
end do

close(3)

open(3, file= "FinalAreas.txt", access="sequential", action = "write", status="replace")
do i = 1, NumSloughs
  write (3,20)  current_area(i)
end do
close(3)

open(3, file= "WetlandDrainageAreas.txt", access="sequential", action = "write", status="replace")
do i = 1, NumSloughs
  write (3,20)  wetland_drainage(i)
end do

close(3)

20 format(1x,f12.3)
40 format(1x, i10,1x,f11.1,1x,f11.1,1x,f11.1,1x,a30,f8.1,1x,f9.1,1x,f9.1,1x,f12.1,1x,f12.1,1x,f8.6,1x,f8.6,1x,f8.6,1x,f8.6)

end program event

