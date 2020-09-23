c     ##########################################################
c     # Program for vector limiting for cutcells               #
c     #                                                        #
c     ##########################################################



      subroutine vector_lim_allineq(num_neighb,neighb_x,neighb_y,
     &    neighb_u,center_x,center_y,center_u,output_debug,row_A,
     & Phi_x,Phi_y,LP_Dx,LP_Dy,LP_success,LS_Dx,LS_Dy,bxpt,bypt,
     & do_pos_constr,my_iters)




c :::::::::::::::::::::::::::: VECTOR_LIM ::::::::::::::::::::::
c
c    This routine does:
c     * set up the 2d LP 
c             - build the matrix A and rhs b
c             - initialize x0 and c appropriately
c             - call the subroutine 'Simplex_allineq' -- that's Simplex in allineq form
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     uses the restrictive constraints (don't overshoot and correct direction)
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  INPUT:
c     * num_neighb: number of neighboring cells used for limiting
c     * neighb_x(num_neighb), neighb_y(num_neighb): x- and y-coordinates of neighbors
c     * neighb_u(num_neighb): cell value of neighbors
c     * center_x,center_y,center_u: data for central cell
c     * output_debug: how much output do you want?
c     * row_A,column_A: dimension of constraint matrix A in LP
c     * LS_Dx,LS_Dy: least square gradient
c     * bxpt,bypt: x- & y-coord of the boundary edge midpoint
c     * do_pos_constr: if true, then add pos constraint for boundary edge midpoint
c  
c  OUTPUT:
c     * Phi_x,Phi_y: scalar factors that limit LS_Dx to LP_Dx (same for y)
c     * LP_Dx,LP_Dy: results of LP: limited slopes
c     * LP_success: was the Simplex algorithm successful?
c     * my_iters: number of iterations needed in the Simplex-algorithm
c
c :::::::::::::::::::::::::::::::::::::::;:::::::::::::::::::::::



c     ::::::::::::::::::::::::::::
c     Declare variables
c     ::::::::::::::::::::::::::::

      implicit none  ! force yourself to properly declare all var

c     Input
      integer num_neighb
      double precision neighb_x(num_neighb), neighb_y(num_neighb) 
      double precision neighb_u(num_neighb)
      double precision center_x, center_y, center_u
      logical output_debug
      integer row_A ! number of constraints in LP, ie number of rows in matrix
      ! number of rows: 2
      double precision LS_Dx, LS_Dy, bxpt, bypt
      logical do_pos_constr

c     Output
      double precision Phi_x, Phi_y 
      double precision LP_Dx, LP_Dy 
      logical LP_success
      integer my_iters

c     Local variables: check num_neighb not too big
      integer max_num_neighb
c      parameter (max_num_neighb = 5)
      parameter (max_num_neighb = 8)
      logical output_debug_here
      logical temp_output

c     Local variables: for least square problem
c      double precision S(num_neighb,2) ! matrix for LS problem
c      double precision S_rhs(num_neighb),rhs_LS(2) ! rhs's for LS problem
c      double precision LS_Dx, LS_Dy    ! solution of LS problem
c      double precision normal_matrix(2,2) ! S'*S
c      character   TRANSA, TRANSB
c      integer ind_perm(2)  ! permutation array coming back from solving matrix
c      integer INFO ! was subroutine of LAPACK successful?

c     Local variables: Formulating the 2d LP
      double precision A(row_A,2)   ! matrix of constraints for LP
      double precision b(row_A)  ! rhs for LP
      double precision x0(2)  ! number of unknowns
      integer W(2) ! working set
      double precision c(2) ! vector c in objective function
      double precision x(2)  ! solution vector of LP

c     Local variables: Loop variables
      integer i,j



c     ::::::::::::::::::::::::::::
c     Declare finished
c     ::::::::::::::::::::::::::::


 902     format (A,'\t',2(E25.15))
 901     format (A,'\t',3e25.15)
 903     format (A,'\t',3(e25.15,'\t'))
 904     format (A,'\t',4(e25.15,'\t'))


c     Check whether num_neighb < max_num_neighb
      if (num_neighb .GT. max_num_neighb) then
         write(*,*) "\n##########\n"
         write(*,901) "num_neighb too big: num_neighb is\t", num_neighb
         write(*,*) "\n##########\n"
         return
      endif


      if (output_debug) then
         write(*,*) "############\n Given data \n ############"
         if (num_neighb .EQ. 3) then
            write(*,903) "neighb_x", neighb_x
            write(*,903) "neighb_y", neighb_y
            write(*,903) "neighb_u", neighb_u
         elseif (num_neighb .EQ. 4) then
            write(*,904) "neighb_x", neighb_x
            write(*,904) "neighb_y", neighb_y
            write(*,904) "neighb_u", neighb_u
         endif
         write(*,902) "center_x, center_y", center_x, center_y
         write(*,901) "center_u", center_u
      endif


c     #####################################
c     Set up and solve least square problem
c     #####################################


cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     LS solution is passed in
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc


cc     IN future: try to do that by hand
c      
c      do 20 i = 1, num_neighb
c         S(i,1) = neighb_x(i) - center_x
c         S(i,2) = neighb_y(i) - center_y
c  20  continue
c      do 30 i = 1, num_neighb
c         S_rhs(i) = neighb_u(i) - center_u
c  30  continue
c      if (output_debug) then
c         write(*,*) "############\n LS problem \n ############"
c         do 31 i = 1, num_neighb
c            write(*,902) "S", (S(i,j),j=1,2)
c 31      continue
c         if (num_neighb .EQ. 3) then
c            write(*,903) "RHS of LS problem", S_rhs
c         elseif (num_neighb .EQ. 4) then
c            write(*,904) "RHS of LS problem", S_rhs
c         endif
c      endif
c
cc     calculate S'*S
c      TRANSA = 'T'
c      TRANSB = 'N'
c      call DGEMM ( TRANSA, TRANSB, 2, 2, num_neighb, 1.0d0, S, 
c     &        num_neighb,S, num_neighb, 0.0d0, normal_matrix,2 )
c
c      if (output_debug) then
c          do 32 i = 1, 2
c            write(*,902) "S'*S", (normal_matrix(i,j),j=1,2)
c 32      continue
c      endif
c
cc     calculate S'*r
c      rhs_LS(1) = 0.0d0
c      rhs_LS(2) = 0.0d0
c      do 33 i = 1, num_neighb
c         rhs_LS(1) = rhs_LS(1) + S(i,1)*S_rhs(i)
c         rhs_LS(2) = rhs_LS(2) + S(i,2)*S_rhs(i)
c 33   continue
c      if (output_debug) then
c         write(*,902) "S'*RHS", rhs_LS
c      endif
c
cc     dgesv: lapack routine for solving Ax=b for x
cc            On exit, if INFO = 0, then (LS_Dx,LS_Dy) = rhs_LS
c      call DGESV( 2,1,normal_matrix,2,ind_perm,rhs_LS,2, INFO)
cc     NOTE: normal_matrix and rhs_LS are lost after this call
c
c      if (INFO.EQ.0) then
c         LS_Dx = rhs_LS(1)
c         LS_Dy = rhs_LS(2)
c      else
c         write(*,*) "\n Error in subroutine dgesv ", 
c     &      "when solving LS problem !!!!!!!!!!!!\n"
c         write(*,*) "Info comes back to", INFO
c         stop
c      endif
c      if (output_debug) then
c         write(*,902) "Solution of LS problem", LS_Dx, LS_Dy
c      endif



c     #######################
c     Formulate 2d LP
c     #######################

         if (output_debug) then
            write(*,*) "############\n Formulate 2d LP \n ############"
         endif

c            ::::::::::::::::::::::
c            initialize A and b
c            ::::::::::::::::::::::

c            Initialize with 0
             do 60 j = 1, 2
                do 70 i = 1, row_A
                   A(i,j) = 0.0d0
  70            continue
  60         continue
             do 80 i = 1, row_A
                b(i) = 0.0d0
  80         continue
c            Add identity matrix in A and the 1's for <=1 condition
             A(1,1) = 1.0d0
             A(2,2) = 1.0d0
             A(3,1) = -1.0d0
             A(4,2) = -1.0d0
             b(3) = -1.0d0
             b(4) = -1.0d0

c            Add data constraints
             do 100 i = 1, num_neighb
                if ((neighb_u(i)-center_u).GT.0) then
                   A(2*i+3,1) = -LS_Dx*(neighb_x(i)-center_x)
                   A(2*i+3,2) = -LS_Dy*(neighb_y(i)-center_y)
                   A(2*i+4,1) = LS_Dx*(neighb_x(i)-center_x)
                   A(2*i+4,2) = LS_Dy*(neighb_y(i)-center_y)
                   b(2*i+3) = center_u-neighb_u(i) 
                   b(2*i+4) = 0.0d0 
               else
                   A(2*i+3,1) = LS_Dx*(neighb_x(i)-center_x)
                   A(2*i+3,2) = LS_Dy*(neighb_y(i)-center_y)
                   A(2*i+4,1) = -LS_Dx*(neighb_x(i)-center_x)
                   A(2*i+4,2) = -LS_Dy*(neighb_y(i)-center_y)
                   b(2*i+3) = neighb_u(i) - center_u
                   b(2*i+4) = 0.0d0
               endif
  100        continue

c          positivity constraint:
            if (do_pos_constr) then
                  A(row_A,1) = LS_Dx*(bxpt-center_x)
                  A(row_A,2) = LS_Dy*(bypt-center_y)
                  b(row_A) = - center_u
            endif

            x0(1) = 0.0d0
            x0(2) = 0.0d0
            c(1) = -abs(LS_Dx)
            c(2) = -abs(LS_Dy)
            x(1) = 0.0d0
            x(2) = 0.0d0
            W(1) = 1
            W(2) = 2

            temp_output = output_debug   ! more output asked for
            LP_success = .FALSE.

c          call 'Simplex_allineq' solver

           if (output_debug) then
            write(*,*) "\n########\n Call Simplex allineq\n#######\n"
           endif

           my_iters = 50  ! should be set to the correct value in Simplex_allineq if that routine was successful

            call Simplex_allineq(row_A,A, b, c,
     &        x0,W,temp_output,x,LP_success,my_iters)
c            write(*,*) "LP_success",LP_success

            if (LP_success) then
               LP_Dx = LS_Dx*x(1)
               LP_Dy = LS_Dy*x(2)
               Phi_x = x(1)
               Phi_y = x(2)
c               write(*,902) "limited slopes", LP_Dx, LP_Dy
c               write(*,*) "\n ###################### \n"
            endif













      return
      end
    

      
