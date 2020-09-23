c
c -------------------------------------------------------------
c
        subroutine finit(val,nvar,maxip1,maxjp1,
     .               corn1,corn2,hx,hy,irr,mitot,mjtot,lstgrd)
c
       implicit double precision (a-h,o-z)
       dimension val(maxip1,maxjp1,nvar),irr(mitot,mjtot)
       dimension state(4)

       common /order2/ ssw, temp2, temp3
       common /userdt/ cfl,gamma,gamma1,xprob,yprob,alpha,Re,iprob,
     .                 ismp,gradThreshold
       common /skewll/ y0skew,dely
      include "cirr.i"
      include "quadrature.i"
      parameter  (maxgr = 192, maxlv=12)
      common /nodal/ rnode(12,maxgr),node(17,maxgr),lstart(maxlv),
     *  newstl(maxlv),
     *  listsp(maxlv),intrat(maxlv),tol,bzonex,bzoney,mstart,ndfree,
     *  lfine,iorder,mxnest,kcheck,lwidth,
     *  maxcl, graf, lhead


      if(ssw .eq. -10 .or. ssw .eq. -1 .or. ssw .eq. 1) then
       ! pointwise
       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 199 i = 1, maxip1
       do 199 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          if(kirr .eq. -1) cycle

          if (kirr .ne. -1 .and. kirr .ne. lstgrd) then
             xcen = xcirr(kirr)
             ycen = ycirr(kirr)
          else if (kirr .eq. lstgrd) then
             xcen = corn1 + (dfloat(i)-1.5d0)*hx
             ycen = corn2 + (dfloat(j)-1.5d0)*hy
          endif

          call f(state, xcen, ycen, 0.d0, iprob)
          val(i,j,:) = state(:)
 199      continue
          return
      endif









       time = 0.d0
       pi = 3.14159265358979d0



       xlow = corn1 - lwidth*hx
       ylow = corn2 - lwidth*hy
       do 100 i = 1, maxip1
       do 100 j = 1, maxjp1
          iadj = i + lwidth - 1
          jadj = j + lwidth - 1
          kirr = irr(iadj,jadj)
          val(i,j,:) = 1.d0


!
!          if(iadj .eq. 47 .and. jadj .eq. 35) then
!          print *,"here"
!          endif






          if (kirr .eq. lstgrd) then ! only full cells for now
             goto 1
          elseif (kirr .ne. -1 .and. kirr .ne. lstgrd
     .           .and. ihob .eq. 0) then ! cut cells
             goto 2
          elseif (kirr .ne. -1 .and. kirr .ne. lstgrd
     .           .and. ihob .eq. 1) then ! cut cells
             goto 3
          else
            cycle
          endif

1         continue ! whole cell projection
          val(i,j,:) = 0.d0



          xcen1 = corn1 + (dfloat(i)-2.d0)*hx
          ycen1 = corn2 + (dfloat(j)-2.d0)*hy

          xcen3 = corn1 + (dfloat(i)-1.d0)*hx
          ycen3 = corn2 + (dfloat(j)-1.d0)*hy


        do 11 iq = 1, nquadquad

      xval = xcen1 * (1.d0 - rquad(iq))/2.d0
     .     + xcen3 * (1.d0 + rquad(iq))/2.d0
      yval = ycen1 * (1.d0 - squad(iq))/2.d0
     .     + ycen3 * (1.d0 + squad(iq))/2.d0


          call f(state, xval, yval, time, iprob)


          val(i,j,1) = val(i,j,1) + wquad(iq) * state(1)
          val(i,j,2) = val(i,j,2) + wquad(iq) * state(2)
          val(i,j,3) = val(i,j,3) + wquad(iq) * state(3)
          val(i,j,4) = val(i,j,4) + wquad(iq) * state(4)


 11   continue
      cycle







 2    continue ! cut cell projection now
          val(i,j,:) = 0.d0
          arr = ar(kirr)
          ivert = 1
          do 20 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  20     continue




          itri = ivert - 3
          idx1 = 1

          tempar = 0.d0
          do 21 it = 1, itri ! for each  triangle
            idx2 = it + 1
            idx3 = it + 2

            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)




            do 22 itq = 1,ntriquad

                xval = x1 * rtri(itq) + x2 * stri(itq)
     .              +  x3 * (1.d0-rtri(itq)-stri(itq))
                yval = y1 * rtri(itq) + y2 * stri(itq)
     .              +  y3 * (1.d0-rtri(itq)-stri(itq))


                  call f(state, xval, yval, time, iprob)


              val(i,j,1) = val(i,j,1) + (artri/arr)*wtri(itq) * state(1)
              val(i,j,2) = val(i,j,2) + (artri/arr)*wtri(itq) * state(2)
              val(i,j,3) = val(i,j,3) + (artri/arr)*wtri(itq) * state(3)
              val(i,j,4) = val(i,j,4) + (artri/arr)*wtri(itq) * state(4)
  22        continue ! for each quadrature point on each triangle

!            print *,val(i,j,1)
  21      continue ! for each triangle
      cycle

 3    continue ! high order cut cell projection now
          val(i,j,:) = 0.d0
          arr = ar_ho(kirr)
          ivert = 1
          do 210 while (poly(ivert+1,1,kirr) .ne. -11.)
            ivert = ivert + 1
  210     continue

          itri = ivert - 3
          idx1 = 1
          tempar = 0.d0
          do 111 it = 1, itri-1 ! for each  triangle except the final one
            idx2 = it + 1
            idx3 = it + 2
            x1 = poly(idx1,1,kirr)
            y1 = poly(idx1,2,kirr)

            x2 = poly(idx2,1,kirr)
            y2 = poly(idx2,2,kirr)

            x3 = poly(idx3,1,kirr)
            y3 = poly(idx3,2,kirr)

            artri = area(x1, x2, x3, y1, y2, y3)
            do 112 itq = 1,ntriquad
             xval = x1 * rtri(itq) + x2 * stri(itq)
     .           +  x3 * (1.d0-rtri(itq)-stri(itq))
             yval = y1 * rtri(itq) + y2 * stri(itq)
     .           +  y3 * (1.d0-rtri(itq)-stri(itq))
             call f(state, xval, yval, time, iprob)
             val(i,j,:) = val(i,j,:) + (artri/arr)*wtri(itq) * state(:)
  112        continue ! for each quadrature point on each triangle
  111      continue ! for each triangle



!      du = (artri/arr)
      if(ssw .eq. 2  .or. ssw .eq. -2) then
      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(3,1,kirr)
      y3 = bdry(3,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)

      do nq = 1,ntriquad_ho
      da = dJ2(x1,x2,x3,x4,y1,y2,y3,y4,nq)
      xc = rs2xy_2(x1,x2,x3,x4,nq)
      yc = rs2xy_2(y1,y2,y3,y4,nq)
      call f(state, xc, yc, time, iprob)
      val(i,j,:) = val(i,j,:) + (da/arr)*wtri_ho(nq) * state(:)
!      du = du + (da/arr)*wtri_ho(nq)
      enddo
!      print *,du


      elseif(ssw .eq. 3 .or. ssw .eq. -3) then



      x1 = poly(ivert-2,1,kirr)
      y1 = poly(ivert-2,2,kirr)
      x2 = bdry(1,1,kirr)
      y2 = bdry(1,2,kirr)
      x3 = bdry(4,1,kirr)
      y3 = bdry(4,2,kirr)
      x4 = bdry(2,1,kirr)
      y4 = bdry(2,2,kirr)
      x5 = bdry(3,1,kirr)
      y5 = bdry(3,2,kirr)

!        du = 0.d0
      do nq = 1,ntriquad_ho
      da = dJ3(x1,x2,x3,x4,x5,y1,y2,y3,y4,y5,nq)
      xc = rs2xy_3(x1,x2,x3,x4,x5,nq)
      yc = rs2xy_3(y1,y2,y3,y4,y5,nq)
      call f(state, xc, yc, time, iprob)
      val(i,j,:) = val(i,j,:) + (da/arr)*wtri_ho(nq) * state(:)
!      du = du + (da/arr)*wtri_ho(nq)
      enddo
      endif

!      if(kirr .eq. 55) then
!      print *, "here"
!      endif
      cycle

 100   continue




      end

      function area(x1, x2, x3, y1, y2, y3)
        implicit double precision (a-h,o-z)

      area = dabs(0.5d0 * ( (x2 - x1) * (y3 - y1)
     .                    - (x3 - x1) * (y2 - y1) ))

      return
      end


