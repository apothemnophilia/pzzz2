      program main
      implicit none
      call input
      call table
      end


      subroutine input
      implicit none
      common/x/ xmin, xmax, xstep
      common/y/ ymin, ymax, ystep
      DOUBLE PRECISION xmin, xmax, xstep, ymin, ymax, ystep

      open(1, file='data.txt', status='old', err=10)
      read(1, *, err=20, end=20) xmin, xmax, xstep
      read(1, *, err=20, end=20) ymin, ymax, ystep
      close(1)

      if (xmin .gt. xmax) then
        write(*,*) 'Error: xmin > xmax'
        stop
      end if
      if (ymin .gt. ymax) then
        write(*,*) 'Error: ymin > ymax'
        stop
      end if
      if (xstep .le. 0.0d0) then
        write(*,*) 'Error: xstep <= 0'
        stop
      end if
      if (ystep .le. 0.0d0) then
        write(*,*) 'Error: ystep <= 0'
        stop
      end if

      write(*, *) 'Input file read successfully.'
      return

10    write(*, *) 'Error: cannot open data.txt'
      stop
20    write(*, *) 'Error: invalid format in data.txt'
      stop
      end


      subroutine table
      implicit none
      common/x/ xmin, xmax, xstep
      common/y/ ymin, ymax, ystep
      DOUBLE PRECISION xmin, xmax, xstep, ymin, ymax, ystep

      DOUBLE PRECISION curr_x, curr_y, sum_deg, nearest_pole
      DOUBLE PRECISION rad_sum, sinval, min_step, pi8
      real*4 f_out
      integer*4 nx, ny, i, j, k
      integer*4 bl_size, nbl, col_in_block
      integer*4 uniq_total, uniq_cnt
      character*13 str_x, prev_str_x, str_y, prev_str_y
      logical first_y_done

      pi8     = 3.14159265358979323846d0
      bl_size = 5


      nx = int((xmax - xmin) / xstep) + 1
      if (abs(xmin + dble(nx-1)*xstep - xmax) .gt.
     &    abs(xstep)*1.0d-5) nx = nx + 1

      ny = int((ymax - ymin) / ystep) + 1
      if (abs(ymin + dble(ny-1)*ystep - ymax) .gt.
     &    abs(ystep)*1.0d-5) ny = ny + 1

      min_step = min(abs(xstep), abs(ystep))

      uniq_total = 0
      prev_str_x = ' '
      do i = 1, nx
          curr_x = xmin + dble(i-1) * xstep
          if (i .eq. nx) curr_x = xmax
          if (abs(curr_x) .lt. abs(xstep) * 0.5d0) curr_x = 0.0d0
          write(str_x, '(e11.4)') curr_x
          if (str_x .ne. prev_str_x) then
              uniq_total = uniq_total + 1
              prev_str_x = str_x
          end if
      end do

      nbl = uniq_total / bl_size
      if (mod(uniq_total, bl_size) .ne. 0) nbl = nbl + 1

      open(2, file='table.txt', status='unknown')

      do k = 1, nbl

          col_in_block = min(bl_size, uniq_total - (k-1)*bl_size)

          write(2, 100) k, nbl
100       format('| BLOCK', i4, ' OF', i4, ' |')

          write(2, 101)
101       format('|    y\x    |', $)

          uniq_cnt  = 0
          prev_str_x = ' '
          do i = 1, nx
              curr_x = xmin + dble(i-1) * xstep
              if (i .eq. nx) curr_x = xmax
              if (abs(curr_x) .lt. abs(xstep) * 0.5d0) curr_x = 0.0d0
              write(str_x, '(e11.4)') curr_x
              if (str_x .ne. prev_str_x) then
                  uniq_cnt = uniq_cnt + 1
                  if (uniq_cnt .gt. (k-1)*bl_size .and.
     &                uniq_cnt .le. k*bl_size) then
                      write(2, 102) curr_x
102                   format(e11.4, '|', $)
                  end if
                  prev_str_x = str_x
              end if
          end do
          write(2, *)

          write(2, 103)
103       format('|-----------|', $)
          do i = 1, col_in_block
              write(2, 104)
104           format('-----------|', $)
          end do
          write(2, *)

          first_y_done = .false.
          prev_str_y   = ' '

          do j = 1, ny
              curr_y = ymin + dble(j-1) * ystep
              if (j .eq. ny) curr_y = ymax
              if (abs(curr_y) .lt. abs(ystep) * 0.5d0) curr_y = 0.0d0

              write(str_y, '(e11.4)') curr_y
              if (first_y_done) then
                  if (str_y .eq. prev_str_y) goto 200
              end if
              prev_str_y   = str_y
              first_y_done = .true.

              write(2, 110) curr_y
110           format('|', e11.4, '|', $)

              uniq_cnt  = 0
              prev_str_x = ' '

              do i = 1, nx
                  curr_x = xmin + dble(i-1) * xstep
                  if (i .eq. nx) curr_x = xmax
                  if (abs(curr_x) .lt. abs(xstep)*0.5d0)
     &                curr_x = 0.0d0
                  write(str_x, '(e11.4)') curr_x
                  if (str_x .ne. prev_str_x) then
                      uniq_cnt = uniq_cnt + 1
                      if (uniq_cnt .gt. (k-1)*bl_size .and.
     &                    uniq_cnt .le. k*bl_size) then

                          sum_deg = curr_x + curr_y
                          nearest_pole =
     &                        dnint(sum_deg / 180.0d0) * 180.0d0
                          if (abs(sum_deg - nearest_pole) .lt.
     &                        min_step * 0.5d0)
     &                        sum_deg = nearest_pole

                          rad_sum = sum_deg * pi8 / 180.0d0
                          sinval  = dsin(rad_sum)

                          if (abs(sinval) .lt. 1.0d-12) then
                              write(2, 120)
120                           format('    inf    |', $)
                          else
                              f_out = real(1.0d0 / sinval)
                              if (abs(f_out) .gt. 9.999e37) then
                                  write(2, 120)
                              else
                                  write(2, 121) f_out
121                               format(e11.4, '|', $)
                              end if
                          end if

                      end if
                      prev_str_x = str_x
                  end if
              end do
              write(2, *)

              if (j .lt. ny) then
                  write(2, 103)
                  do i = 1, col_in_block
                      write(2, 104)
                  end do
                  write(2, *)
              end if

200           continue
          end do

          write(2, 103)
          do i = 1, col_in_block
              write(2, 104)
          end do
          write(2, *)

          write(2, 130) k, nbl
130       format('| END BLOCK', i4, ' OF', i4, ' |')
          write(2, *)

      end do

      close(2)
      write(*, *) 'Table created: table.txt'
      end
