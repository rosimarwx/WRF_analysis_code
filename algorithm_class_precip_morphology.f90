subroutine class_precip_morphology(ni,nj,nk,reflectivity,www,varout,DY,DX)
  implicit none
  integer, intent(IN) :: nk, nj, ni
  real, intent(IN)    :: DX, DY
  real, intent(IN)    :: reflectivity(ni,nj,nk), www(ni,nj,nk)
  real, intent(OUT)   :: varout(ni,nj)
  integer :: kk, jj, ii, jpts, ipts, j1, j2, i1, i2, iii, jjj
  real :: Zbg(ni,nj,nk), cnt(ni,nj,nk), convrad, dist, dummyout(ni,nj,nk)

  !OUTPUT: dummyout = 1 if convective, = 2 if stratiform, = 3 other, = 0 if no rain

  dummyout(:,:,:) = 0.
  varout(:,:) = 0.
  Zbg = 0.
  cnt = 0.

  do kk = 1, nk
    do jj = 1, nj; do ii = 1, ni
      !First pass: dBZ > 45 = convective 
      if (reflectivity(ii,jj,kk) .gt. 45.) then
        dummyout(ii,jj,kk) = 1.
      else
        !Second pass: if averaged w < 0.1 m/s then it's unlikely to be convective
        if (sum(www(ii,jj,:)/nk) .lt. 0.5) then
          cycle
        else
          !Second pass: if dBZ < 46, check background reflectivity of nonzero echoes
          jpts = int(11/DY)
          jpts = int(11/DX)
          j1 = max(jj-jpts, 1)
          i1 = max(ii-ipts, 1)
          j2 = min(nj,jj+jpts)
          i2 = min(ni,ii+ipts)
          do jjj = j1, j2; do iii = i1, i2
            dist = sqrt(((ii-iii)*DX)**2.+((jj-jjj)*DY)**2.)
            if (dist .le. 11.) then
              Zbg(ii,jj,kk) = Zbg(ii,jj,kk)+reflectivity(iii,jjj,kk)
              cnt(ii,jj,kk) = cnt(ii,jj,kk)+1.
            end if
          end do; end do

          Zbg(ii,jj,kk) = Zbg(ii,jj,kk)/cnt(ii,jj,kk)
      
          if (Zbg(ii,jj,kk) .lt. 0. .and. (reflectivity(ii,jj,kk) - Zbg(ii,jj,kk)) .gt. 10.) then
            dummyout(ii,jj,kk) = 1.
          else if (Zbg(ii,jj,kk) .ge. 0. .and. Zbg(ii,jj,kk) .lt. 45. &
                  .and. (reflectivity(ii,jj,kk) - Zbg(ii,jj,kk)) .gt. &
                  (15.-Zbg(ii,jj,kk)*Zbg(ii,jj,kk)/135.) ) then
            dummyout(ii,jj,kk) = 1.
          else if (Zbg(ii,jj,kk) .ge. 45. .and. (reflectivity(ii,jj,kk) - Zbg(ii,jj,kk)) .gt. 0.) then
            dummyout(ii,jj,kk) = 1.
          end if     

        end if

      end if
    end do; end do

    !Third pass - for each grid point identified as convective above, all surrounding grid points within an intensity-dependent convective radius are also convective
    do jj = 1, nj; do ii = 1, ni
      if (dummyout(ii,jj,kk) .eq. 1.) then !identified convective
        !define radius based on intensity
        if (Zbg(ii,jj,kk) .le. 25.) then
          convrad = 1.
        else if (Zbg(ii,jj,kk) .gt. 25. .and. Zbg(ii,jj,kk) .le. 30.) then
          convrad = 2.
        else if (Zbg(ii,jj,kk) .gt. 30. .and. Zbg(ii,jj,kk) .le. 35.) then
          convrad = 3.
        else if (Zbg(ii,jj,kk) .gt. 35. .and. Zbg(ii,jj,kk) .le. 40.) then
          convrad = 4.
        else if (Zbg(ii,jj,kk) .gt. 40.) then
          convrad = 5.
        end if

        jpts = int(convrad/DY)
        jpts = int(convrad/DX)
        j1 = max(jj-jpts, 1)
        i1 = max(ii-ipts, 1)
        j2 = min(nj,jj+jpts)
        i2 = min(ni,ii+ipts)
        do jjj = j1, j2; do iii = i1, i2
          dist = sqrt(((ii-iii)*DX)**2.+((jj-jjj)*DY)**2.)
          if (dist .le. convrad) then    
            dummyout(ii,jj,kk) = 1.
          end if    
        end do; end do       
      end if
    end do; end do

  end do

  !Last pass - put everything together. Any point identified as convective in any grid point gets 
  !convective. Otherwise stratiform, no rain, or other
  do jj = 1, nj; do ii = 1, ni; do kk = 1, nk
    if (dummyout(ii,jj,kk) .eq. 1.) then
      varout(ii,jj) = 1.
      exit
    else if (dummyout(ii,jj,kk) .eq. 0. .and. (reflectivity(ii,jj,kk) .gt. 10.) ) then
      varout(ii,jj) = 2.
    else if (dummyout(ii,jj,kk) .eq. 0. .and. (reflectivity(ii,jj,kk) .ge. 0. .and. &
               reflectivity(ii,jj,kk) .le. 10.) ) then
      varout(ii,jj) = 3.
    end if    
  end do; end do; end do

  return 
  end subroutine class_precip_morphology
