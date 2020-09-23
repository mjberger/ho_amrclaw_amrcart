c
c ------------------------------------------------------------------
c
       subroutine evalU(Uout,deltax,deltay,dx,dy,q,qx,qy,
     &                  qxx,qxy,qyy,i,j,kirr,mitot,mjtot,nvar)

       use amr_module

       implicit double precision (a-h,o-z)

       dimension q(nvar,mitot,mjtot)
       dimension qx(nvar,mitot,mjtot), qy(nvar,mitot,mjtot)
       dimension qxx(nvar,mitot,mjtot), qyy(nvar,mitot,mjtot)
       dimension qxy(nvar,mitot,mjtot)
       dimension Uout(nvar)

       ! this assumes ssw = 2 for now, meaning quadratic poly in each cell
       ! and linear boundary segments

       ! deltax,y is distance from centroid to edge
       ! dx,dy is mesh width


       Uout(:) = q(:,i,j) + deltax*qx(:,i,j)+ deltay*qy(:,i,j)
     &          + 0.5d0*qxx(:,i,j)*(deltax**2-(dx**2)*poly(8,1,kirr)  )
     &          +       qxy(:,i,j)*(deltax*deltay-dx*dy*poly(9,1,kirr))
     &          + 0.5d0*qyy(:,i,j)*(deltay**2-(dy**2)*poly(10,1,kirr) )

      return
      end
