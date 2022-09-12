   !
   !                        libxmd
   !                      version 2.1.2
   !
   !                 (c) 2011   M. Ibaraki
   !
   !        version 2.0     April 30, 2011
   !                    include realloc and ialloc with arrays for semi- list.append() function in python
   !                    must be replaced with c program
   !        version 2.1     May 8, 2011
   !                    no scaling is required for red-black ordering (reduced system)
   !        version 2.1.1   May 12, 2011
   !                    minor changes
   !        version 2.1.2   January 23, 2012
   !                    resolved cgstab and orthomin initialization issues
   !
   !
module xmdalloc
contains
   subroutine realloc(a, n)
   !         a: array
   !         n: new size
      implicit none
      double precision, allocatable, dimension(:) :: a, tmp
      integer :: n, nold, ierr
      nold = size(a)
      if (nold > n) return
      allocate (tmp(nold+200000), stat=ierr)
      if (ierr /= 0) stop "allocate error"
      tmp(1:nold) = a(1:nold)
      deallocate (a)
      allocate (a(n+200000), stat=ierr)
      if (ierr /= 0) stop "allocate error"
      a(1:nold) = tmp(1:nold)
   end subroutine realloc

   subroutine ialloc(a, n)
   !         a: array
   !         n: new size
      implicit none
      integer, allocatable, dimension(:) :: a, tmp
      integer :: n, nold, ierr
      nold = size(a)
      if (nold > n) return
      allocate (tmp(nold+200000), stat=ierr)
      if (ierr /= 0) stop "allocate error"
      tmp(1:nold) = a(1:nold)
      deallocate (a)
      allocate (a(n+200000), stat=ierr)
      if (ierr /= 0) stop "allocate error"
      a(1:nold) = tmp(1:nold)
   end subroutine ialloc

end module xmdalloc





module xmdcmn
   !
   !     miunit      unit number for output
   !     miout       switch for information message
   !                  <  0 nothing printed
   !                  >= 0 print all
   !
   integer ::  miunit, miout

end module xmdcmn



module xmdmatrix

   !     icolour()        nn       poiter array for red-black ordered matrix
   !     RBorder()       nn        Red-Black ordering vector
   !     iblackend()       nn      end markder for black nodes in a row
   !     iaf()        nblack+1     YSMP array for factored matrix
   !     jaf()         njaf        YSMP array for factored matrix
   !     idiagf()     nblack       diagonal location for factored matrix
   !     af()          njaf        factored matrix (each row of af contains a row L\U)
   !                               where A = LU
   !     lorder()      nn          ordering vector  lorder(new_number) = old_number
   !     njaf                      size of jaf() and af()
   !     nblack                    number of black nodes

   integer, dimension(:), allocatable :: icolour, RBorder, iblackend,&
   &iaf, idiagf, lorder
   integer :: njaf, nblack

   integer, dimension(:), allocatable :: jaf
   double precision, dimension(:), allocatable :: af

end module xmdmatrix





module xmdSymFacLev
contains
   subroutine xmdsfacl(iaf, jafwk, idiagf, ia, ja,&
   &icolour, RBorder, iblackend,&
   &n, nja, njaf, level, nblack, ierr)

      use xmdcmn
      use xmdalloc

      implicit none
      integer :: ia(n+1), ja(nja), idiagf(nblack), iaf(nblack+1),&
      &icolour(n), RBorder(n), iblackend(n),&
      &n, nja, njaf, ierr, level, nblack

   !     integer, dimension(:), allocatable :: jaf
   !
   !
   !     incomplete LU factorization
   !
   !     jaf() is wk array, on exit first iaf(n+1)-1 entries
   !     are ja() factors
   !
   !
   !     output: iaf, jaf, idiagf, njaf
   !
   !
   !     lorder(new_order)=old_order
   !
   !     brute force factor
   !
   !     local variables
   !
   !
      integer :: itemp, ii, first, iold, iend, num, next,&
      &iblck, id, kk, nred, ierror
      integer :: maxint = 999999

      integer, dimension(:), allocatable :: levptr1, list

      integer, dimension(:), allocatable :: jafwk, levptr, jafwk0
      if (allocated(jafwk)) deallocate (jafwk)
      allocate( jafwk(1), levptr(1), stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdsfacl) =="

      allocate( levptr1(n), list(nblack), stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdsfacl) =="

   !     initialize arrays

      iaf(1:nblack+1) = 0
      list(1:nblack) = 0
      idiagf(1:nblack) = 0
      levptr1(1:n) = maxint

      iaf(1) = 1

      itemp = 0                    ! current position in jaf

      nred = n - nblack

      BlackNodes: do iblck = 1, nblack       ! do it in NEW number!!

         iold = RBorder(iblck)       ! black node number in old numbering

   !     load row i of L/U into list  -- black node

         iend = itemp

   !     diagonal is black

         id = iblck
         iend = iend + 1

         call ialloc(jafwk, iend)

         jafwk(iend) = id

   !       off-diagonals - black node

         do ii = ia(iold)+1, iblackend(iold)
            id = icolour(ja(ii))
            iend = iend + 1
            call ialloc(jafwk, iend)
            jafwk(iend) = id
         enddo

   !     load row i of L/U into list  -- red node elimination

         do ii = iblackend(iold)+1, ia(iold+1)-1 ! off-diagonals and red

            id = ja(ii)

   !       ignore black node or diagonal node

            do kk = ia(id)+1, ia(id+1)-1   ! off-diagonals

               iend = iend + 1
               call ialloc(jafwk, iend)
               jafwk(iend) = icolour(ja(kk))

            enddo

         enddo

         num = iend-itemp

   !       sort entries

         call xmdshell(jafwk(itemp+1), num)

         do ii = itemp+1, iend
            levptr1( jafwk(ii) ) = 0

   !      new fill-in caused by red node elimination -> level = 1

            if (icolour(jafwk(ii)) < 0) levptr1( jafwk(ii) ) = 1
         enddo


         first = jafwk(itemp+1)
         do ii = itemp+1, iend-1
            list(jafwk(ii)) = jafwk(ii+1)   ! link list containing next column number
         enddo
         list(jafwk(iend)) = n+1         ! end marker

         call xmdmrgl(iblck, iaf, jafwk, idiagf, list, size(jafwk),&
         &nblack,&
         &first, level, levptr, levptr1, size(levptr))

         next = first
100      continue
         if (next /= n+1) then
            itemp = itemp+1
            call ialloc(jafwk, itemp)
            jafwk(itemp) = next

            call ialloc(levptr, itemp)
            levptr(itemp) = levptr1( next )

            levptr1( next ) = maxint

            if (next == iblck) idiagf(iblck) = itemp

            next = list(next)   ! link list
            goto 100
         endif
         iaf(iblck+1) = itemp+1
         if (idiagf(iblck) == 0) then
            ierr  =  3
            write (miunit,200) iblck
            return
         endif

      enddo BlackNodes


200   format ('  error in xmdsfacl'/&
      &'    no diagonal in L\U:  row number',i8)

   !     minjaf = iaf(nblack+1)-1   ! size of af and jaf

      njaf = itemp
   !     print *, 'actual size of af, jaf arrays (njaf): ', njaf

   !     the size of jafwk is bigger than njaf...

      allocate (jafwk0(njaf))
      jafwk0(1:njaf) = jafwk(1:njaf)
      deallocate (jafwk)
      allocate (jafwk(njaf))

      jafwk(1:njaf) = jafwk0(1:njaf)

   !     deallocate(levptr, levptr1, list, jafwk0)
   !     nullify(jafwk0, levptr)

   end subroutine xmdsfacl
end module xmdSymFacLev







module xmdSymFacDrop
contains

   subroutine xmdsfacd(a, b, afwk, epsrn, iaf, jafwk, idiagf,&
   &ia, ja, icolour, RBorder, iblackend,&
   &n, nja, njaf, level, nblack, ierr)


      use xmdcmn
      use xmdalloc
      integer n, nja, njaf, ierr, level, nblack
      integer ia(n+1), ja(nja), idiagf(nblack),&
      &iaf(nblack+1),&
      &icolour(n), RBorder(n), iblackend(n)

      double precision a(nja), b(n), epsrn
   !
   !
   !     incomplete LU factorization which uses the combination of
   !     level and drop tolerance
   !
   !
   !     input:
   !       a()            coefficient matrix
   !       b()            RHS vector
   !       epsrn          drop tolerance
   !       ia(), ja()     YSMP pointers for a
   !       icolour()        pointer for node attribute
   !                > 0  -> black node number
   !                < 0  -> red   node number
   !       iblackend()       pointer to indicate the end of black node in a row
   !       RBorder()       permutation vector for black-red
   !          - RBorder(1:nblack)
   !             original_node_number = RBorder(black_node_number)
   !
   !          - RBorder(nblack+1:neq)
   !             original_node_number = RBorder(red_node_number)
   !                                  note: decrease in numbers
   !       n              total number of nodes
   !       nja            size of ja()
   !       njaf           size of af() and jaf()
   !       lsize          size of levptr ( >= njaf+nblack )
   !       level          level of ILU
   !       nblack         number of black nodes
   !     output:
   !       af()           factored coefficient matrix
   !       iaf(), jaf()   YSMP pointers for af
   !       idiagf()       diagonal pointer for af
   !       ierr           error flag
   !     temporary:
   !       row()          real temp work array
   !       list()         integer temp work array
   !       levptr()       integer temp work array and devide this as
   !                      [njaf][n]
   !                       ^^^^  ^
   !                [level in af][level in row i]
   !
   !
   !     jaf() is wk array, on exit first iaf(n+1)-1 entries
   !     are ja() factors
   !
   !     brute force factor
   !
   !     local variables
   !
      integer :: itemp, ii, first, iold, iend, num, next,&
      &maxint, iblck, id, kk, idk, idkk, ierror
      double precision temp, epsmin

      double precision, dimension(:), allocatable :: row
      integer, dimension(:), allocatable :: list, levptr1
      parameter (maxint = 999999)

      integer, dimension(:), allocatable :: jafwk, levptr, jafwk0
      double precision, dimension(:), allocatable :: afwk, afwk0

      allocate( levptr(1) )
      if (allocated(jafwk)) deallocate (jafwk)
      if (allocated(afwk)) deallocate (afwk)
      allocate( jafwk(1),  afwk(1) , stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdsfacd) =="


   !     lsize = njaf+nblack

      allocate( row(nblack), stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdsfacd) =="

      allocate( list(nblack), levptr1(nblack), stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdsfacd) =="

      epsmin = 1.0d-300

   !     initialize arrays

      levptr1(1:nblack) = maxint
      iaf(1:nblack+1) = 0
      list(1:nblack) = 0
      idiagf(nblack) = 0
      row(1:nblack) = 0.0d0

      iaf(nblack+1) = 0
      iaf(1) = 1
      itemp = 0                    ! current position in jaf

      do 4 iblck = 1, nblack       ! do it in NEW number!!

         iold = RBorder(iblck)       ! black node number in old numbering

   !     load row i of L/U into list  -- black node

         iend = itemp

   !     diagonal is black

         id = iblck
         iend = iend + 1
         call ialloc(jafwk, iend)
         jafwk(iend) = id

         row(id) = a( ia(iold) )

   !       off-diagonals - black node

         do ii = ia(iold)+1, iblackend(iold)
            id = icolour(ja(ii))
            iend = iend + 1
            call ialloc(jafwk, iend)
            jafwk(iend) = id
            row(id) = row(id) + a(ii)
         enddo

   !     load row i of L/U into list  -- red node elimination

         do ii = iblackend(iold)+1, ia(iold+1)-1 ! off-diagonals and red

            id = ja(ii)

   !       ignore black node or diagonal node

            do kk = ia(id)+1, ia(id+1)-1   ! off-diagonals

               iend = iend + 1
               call ialloc(jafwk, iend)
               jafwk(iend) = icolour(ja(kk))

               idk = ja(kk)
               idkk = icolour(idk)

               row(idkk) = row(idkk) - a(ii) * a(kk) / a( ia(id) )

            enddo

   !      modify {b}

            b(iold) = b(iold) - a(ii)*b(id) / a( ia(id) )

         enddo

         num = iend-itemp

   !       sort entries

         call xmdshell(jafwk(itemp+1), num)

         do ii = itemp+1, iend

            levptr1(jafwk(ii)) = 0

   !      new fill-in caused by red node elimination -> level = 1

            if (icolour(jafwk(ii)) < 0) levptr1(jafwk(ii)) = 1

         enddo


         first = jafwk(itemp+1)
         do ii = itemp+1, iend-1
            list(jafwk(ii)) = jafwk(ii+1)   ! link list containing next column number
         enddo
         list(jafwk(iend)) = n+1         ! end marker


         call xmdmrgd(a, afwk, row, epsrn, iblck, nblack, ia, iaf,&
         &jafwk, idiagf, list, RBorder, nja, size(jafwk),&
         &n, first, level, levptr, levptr1, size(levptr))

         next = first
100      continue
         if (next /= n+1) then
            itemp = itemp+1
            call ialloc(jafwk, itemp)
            jafwk(itemp) = next

            call ialloc(levptr, itemp)
            levptr(itemp) = levptr1(next)
            levptr1(next) = maxint

            if (next == iblck) idiagf(iblck) = itemp

            next = list(next)
            goto 100
         endif
         iaf(iblck+1) = itemp+1
         if (idiagf(iblck) == 0) then
            ierr  =  3
            write (miunit,200) iblck
            return
         endif

   !     inverse of diagonal to avoid small pivots

         temp = row(iblck)
         ii = idiagf(iblck)

         call realloc(afwk, ii)
         afwk(ii) = 1.0d0/(temp + epsmin)

   !
   !     gathering full length array
   !
         do ii = iaf(iblck), iaf(iblck+1)-1

            call realloc(afwk, ii)
            afwk(ii) = row(jafwk(ii))

            row(jafwk(ii)) = 0.0d0
            list(jafwk(ii)) = 0
         enddo

4     continue

200   format ('  error in xmdsfacd'/&
      &'    no diagonal in L\U:  row number',i8)

   !     minjaf = iaf(nblack+1)-1   ! size of af and jaf

      njaf = itemp

   !     print *, 'actual size of af, jaf arrays (njaf): ', njaf

   !     the size of jafwk, afwk is bigger than njaf...

   !     jafwk0 => jafwk   ! ....wk0 are tmp arrays
   !     afwk0 => afwk

      allocate (jafwk0(njaf), afwk0(njaf))
      jafwk0(1:njaf) = jafwk(1:njaf)
      afwk0(1:njaf) = afwk(1:njaf)

      deallocate (jafwk, afwk)
      allocate (jafwk(njaf), afwk(njaf))
      jafwk(1:njaf) = jafwk0(1:njaf)
      afwk(1:njaf) = afwk0(1:njaf)

      deallocate(row, list, levptr, levptr1, jafwk0, afwk0)
   !     nullify(jafwk0, afwk0, levptr)

   end subroutine xmdsfacd

end module xmdSymFacDrop









subroutine xmdordng(ia, ja, lorder, neq, nja, norder, ierr)

   use xmdcmn
   !
   !     compute ordering
   !
   !     input:
   !        ia(), ja()    YSMP arrays
   !        nja           dim of ja()
   !        neq           number of nodes
   !        norder        ordering scheme
   !                      = 0  natural ordering
   !                      = 1  RCM ordering
   !                      = 2  minimum degree
   !     output:
   !        lorder()      ordering vector
   !                         lorder(new_number) = old_number
   !        ierr          error flag
   !
   integer :: neq, nja
   integer :: ia(neq+1), ja(nja), lorder(neq), norder, ierr
   !
   !     local variables
   !
   integer :: i, nsp, iflag, ierror
   integer, dimension(:), allocatable :: iwork0, iwork1

   allocate( iwork0(neq), iwork1(neq+1), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdordng) =="

   ierr = 0

   !  ---------------
   !     get lorder
   !  ---------------
   if (norder == 0) then
   !
   !     natural ordering
   !
      do i = 1, neq
         lorder(i) = i
      enddo

   elseif (norder == 1) then
   !
   !     RCM ordering
   !
      call genrcm(neq, nja, ia, ja,&
      &lorder, iwork0, iwork1)

   !       call genrcm(neq, nja, ia, ja, lorder)

   elseif (norder == 2) then
      deallocate( iwork1 )
      nsp = 3*neq + 4*nja
      allocate( iwork1(nsp), stat = ierror )
      if (ierror /= 0) stop "== not enough memory (xmdordng) =="
   !
   !     use odrv code to get minimum degree ordering
   !
      call odrv(ia, ja, lorder, iwork0, iwork1,&
      &neq, nja, nsp, iflag)

      if (iflag /= 0) then
         write (miunit, 100) iflag
         ierr = 3
         return
      endif
   endif

100 format ('  error in min. degree ordering'/&
   &'    error flag',i8)

   !     deallocate( iwork0, iwork1 )
end subroutine xmdordng




subroutine xmdRedBlack(ia, ja, lorder, icolour, RBorder,&
&iblackend, neq, nja, nblack, ierr, redcdsys)
   !
   !     red-black ordering setup - compute red-black ordering vectors
   !
   !     input:
   !           neq      number of equations
   !           ia()     YSMP array
   !           lorder() ordering vector
   !                      lorder(new_number) = old_number
   !           nja      size of ja array
   !           redcdsys = .true.      reduced system
   !                    = .false.     full system
   !     output:
   !           nblack   number of black node
   !           icolour()  > 0  -> black node number
   !                    < 0  -> red   node number
   !           RBorder()
   !              - RBorder(1:nblack)
   !                 original_node_number = RBorder(black_node_number)
   !
   !              - RBorder(nblack+1:neq)
   !                 original_node_number = RBorder(red_node_number)
   !
   !           iblackend()   end marker for black nodes in a row
   !           ierr     error flag
   !     in-output:
   !           ja()     YSMP array  - order of entries in a row is changed
   !
   !
   !     temp:
   !           ilist()
   !
   use xmdcmn

   implicit none
   integer :: ia(neq+1), ja(nja), lorder(neq), nblack, neq, nja, ierr
   integer :: icolour(neq), RBorder(neq), iblackend(neq)
   logical :: redcdsys
   !
   !     locall variables
   !
   integer :: ibf, nred, ii, i, iHalf, iBlack,iRed, ierror, k
   integer, allocatable, dimension(:) :: ilist
   allocate( ilist(neq) , stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdRedBlack) =="

   ierr = 0
   nred = 0
   nblack = 0
   !
   !  determine red or black nodes accoring to their connections
   !
   !     ilist(ibk) = 0 : undetermined
   !                = 1 : red node
   !                = -1 : black node - neighbour of a red node
   !
   if (.not.redcdsys) then    ! full system
      ilist(1:neq) = -1           ! all nodes are black
      do i  = 1, neq
         ibf = lorder(i)
         call xmdBlackCounter(icolour, RBorder, neq, nblack, ibf)
      enddo
   else                  ! reduced system
      ilist(1:neq) = 0            ! set all nodes are undetermined
   !
   !  check all nodes to find red nodes (surrounded by black nodes)
   !
      do i  = 1, neq

         ibf = lorder(i)

         if (ilist(ibf) == -1) then  ! ibf is a black node
            call xmdBlackCounter(icolour, RBorder, neq, nblack, ibf)
            cycle
         endif

   !  check ibf's neighbour,

         do ii = ia(ibf), ia(ibf+1)-1
            if (ilist( ja(ii) ) == 1) then ! found a red node in neighbour. ibf should be a black node
               call xmdBlackCounter(icolour, RBorder, neq, nblack, ibf)
               cycle
            endif
         enddo

   !  this is a red node!

         ilist(ibf) = 1

   !    neighbouring nodes should be black

         do ii = ia(ibf), ia(ibf+1)-1
            if (ja(ii)  /=  ibf) ilist(ja(ii)) = -1
         enddo

   !    count red node

         nred = nred + 1
         icolour(ibf) = -nred  ! negative # for red node
         RBorder( neq +1 -nred) = ibf  ! store red node info into RBorder

      enddo

   endif


   if (nblack + nred  /=  neq) then
      ierr = 3
      write (miunit,100)
      return
   endif

100 format ('  error in subroutine xmdRedBlack'/&
   &'   sum of n_black and n_red <> n')
200 format ('  need more space in integer temp. array')
300 format ('  error in red-black ordering'/&
   &'    total number of nodes',i8/&
   &'    number of red node',i8/&
   &'    number of black node',i8)






   !
   !     rearrange column numbers
   !
   iHalf = neq/2
   do i = 1, neq

      iBlack = 0
      iRed = 0
      do ii = ia(i)+1, ia(i+1)-1  ! check neighbours

         ibf = ja(ii)

         if (icolour( ibf ) > 0) then        ! black node
            iBlack = iBlack + 1
            if (iBlack > iHalf) then
               ierr = 2
               write (miunit,200)
               return
            endif
            ilist(iBlack) = ibf
         else
            iRed = iRed + 1
            if (iRed > iHalf) then
               ierr = 2
               write (miunit,200)
               return
            endif
            ilist(iHalf + iRed) = ibf
         endif
      enddo


   !     check number of entries

      k = ia(i+1) - ia(i) - 1
      if (k /= iRed+iBlack) then
         ierr = 3
         write (miunit,300) k, iRed, iBlack
         return
      endif

   !     sort entries

      call xmdshell(ilist(1), iBlack)
      call xmdshell(ilist(iHalf+1), iRed)

      k = ia(i)
      do ii = 1, iBlack
         k = k + 1
         ja(k) = ilist(ii)
      enddo

      iblackend(i) = k  ! end marker for black nodes

      do ii = 1, iRed
         k = k + 1
         ja(k) = ilist(iHalf+ii)
      enddo
   !
   enddo

   !     deallocate( ilist )
end subroutine xmdRedBlack







subroutine xmdBlackCounter(icolour, RBorder, neq, nblack, ibf)
   !
   !     count black node
   !
   implicit none
   integer :: icolour(neq), RBorder(neq), neq, nblack, ibf

   nblack = nblack+1
   icolour(ibf) = nblack   ! positive # for black node
   RBorder(nblack) = ibf

end subroutine xmdBlackCounter



subroutine xmdshell(x,n)
   !
   !        shell sort
   !
   implicit none
   integer n, x(n)
   !
   !
   !       return x(i), i=1,..n
   !         sorted in increasing value
   !
   !       local vars
   !
   integer m, maxi, k, j, itemp, xert
   !
   m = n
10 continue
   m = m/2
   if( m  ==  0) return
   maxi = n - m
   do 20 j=1, maxi
      do k = j, 1, -m
         xert = x(k+m)
         if( xert  >=  x(k) ) go to 20
         itemp = x(k+m)
         x(k+m) = x(k)
         x(k) = itemp
      enddo
20 continue
   goto 10
end subroutine xmdshell






subroutine xmdcgstb(a, b, x, af, soln,&
&ctol, rrctol,&
&ia, ja, iaf, jaf, idiagf, RBorder, nblack,&
&nred, n, nja, njaf, nitmax, north, ierr)


   implicit none
   double precision :: a(nja), b(n), af(njaf), x(n), soln(nblack),&
   &ctol, rrctol


   integer :: ia(n+1), ja(nja), iaf(nblack+1), jaf(njaf),&
   &idiagf(nblack), RBorder(n), n, nja, njaf, nitmax, nblack,&
   &nred, north, ierr
   !
   !     local variables
   !
   !     work arrays
   !        q()     nblack
   !        qb()    nblack
   !        aqb()   nblack
   !        reso()  nblack
   !        res()   nblack
   !        t()     nblack
   !        s()     nblack
   !        sb()    nblack
   !       total     9*nblack

   integer :: i, iter, irpnt, ierror
   double precision :: res0, omega, omegah, beta, betah, alpha,&
   &resmax, temp1, temp2, temp, verysmall
   logical :: conv, rescal

   double precision, dimension(:), allocatable :: q, qb, aqb, reso,&
   &res, t, s, sb

   allocate( q(nblack), qb(nblack), aqb(nblack), reso(nblack),&
   &res(nblack), t(nblack), s(nblack),&
   &sb(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdcgstb) =="

   verysmall = 1.0d-300

   q(1:nblack) = 0.0d0
   aqb(1:nblack) = 0.0d0

   irpnt = nblack+1             ! pointer for REDorder
   if (nred < 1) irpnt = 1

   !     get [a]{x} -- note: use t as tmp

   call xmdrbmtv(a, x, t, ia, ja, RBorder(1), RBorder(irpnt),&
   &n, nja, nblack, nred)

   res0 = 0.0d0
   do i = 1, nblack
      temp = b( RBorder(i) ) - t(i)
      t(i) = 0.0d0
      res(i) = temp
      reso(i) = res(i)
      res0 = res0 + temp * temp
   enddo

   res0 = dsqrt(res0)
   omegah = 1.0d0
   beta = 1.0d0
   alpha = 1.0d0

   !     save x into soln and use x as temporary(n) array

   do i = 1, nblack
      soln(i) = x( RBorder(i) )
   enddo

   !      soln(1:nblack) = x( RBorder(1:nblack) )  !rgn

   !
   !     iteration loop
   !
   do iter = 0, nitmax

      betah = 0.0d+0
      do i = 1, nblack
         betah = betah + reso(i)*res(i)
      enddo

   !       omega = (betah/beta)*(omegah/alpha)
      omega = betah/(beta+sign(verysmall, beta)) *&
      &omegah/(alpha+sign(verysmall, alpha))
      beta = betah
      do i = 1, nblack
         q(i) = res(i)+omega*(q(i)-alpha*aqb(i))
      enddo

      call xmdilusl(qb, q, af, iaf, jaf, idiagf, njaf,&
      &nblack)
   !
   !     CGSTAB acceleration
   !
      do i = 1, nblack
         x(RBorder(i)) = qb(i)
      enddo
      call xmdrbmtv(a, x, aqb, ia, ja, RBorder(1), RBorder(irpnt),&
      &n, nja, nblack, nred)

      temp1 = 0.0d+0
      do i = 1, nblack
         temp1 = temp1 + reso(i)*aqb(i)
      enddo

   !       omegah = betah/temp1
      omegah = betah/(temp1+sign(verysmall, temp1))

      do i = 1, nblack
         s(i) = res(i) - omegah*aqb(i)
      enddo

      call xmdilusl(sb, s, af, iaf, jaf, idiagf, njaf, nblack)

      do i = 1, nblack
         x(RBorder(i)) = sb(i)
      enddo
      call xmdrbmtv(a, x, t, ia, ja, RBorder(1), RBorder(irpnt),&
      &n, nja, nblack, nred)

      temp1 = 0.0d+0
      temp2 = 0.0d+0
      do i = 1, nblack
         temp1 = temp1 + t(i)*t(i)
         temp2 = temp2 + t(i)*s(i)
      enddo

   !       alpha = temp2/temp1
      alpha = temp2/(temp1+sign(verysmall, temp1))

      resmax = 0.0d+0
      temp1 = 0.0d+0
      conv = .false.

      do i = 1, nblack
         soln(i) = soln(i) + omegah*qb(i) + alpha*sb(i)
         res(i) = s(i) - alpha*t(i)
         temp1 = dmax1( dabs(omegah*qb(i)), dabs(qb(i)),&
         &dabs(alpha*sb(i)), dabs(sb(i)),&
         &dabs(omegah*qb(i) + alpha*sb(i)), temp1 )
         resmax = resmax + res(i)*res(i)
      enddo

      resmax = dsqrt(resmax)

      conv = temp1 < ctol
      conv = conv.or.(resmax < rrctol*res0)

      if (conv) then
         nitmax = iter+1
   !         deallocate( q, qb, aqb, reso, res, t, s, sb )
         return
      endif
   !
   !     explicit calculation of residual at every ``north'' step
   !
      rescal = mod(iter+1,north) == 0

      if (rescal) then
         do i = 1, nblack
            x(RBorder(i)) = soln(i)
         enddo
         call xmdrbmtv(a, x, qb, ia, ja, RBorder(1), RBorder(irpnt),&
         &n, nja, nblack, nred)
         do i = 1, nblack
            res(i) = b( RBorder(i) ) - qb(i)
         enddo
      endif

   enddo

   !     no convergence

   ierr = -1
   !     deallocate( q, qb, aqb, reso, res, t, s, sb )

end subroutine xmdcgstb



subroutine xmdcnjgd(a, b, x, af, soln,&
&ctol, rrctol,&
&ia, ja, iaf, jaf, idiagf, RBorder, nblack,&
&nred, n, nja, njaf, nitmax, ierr)


   implicit none
   double precision :: a(nja), b(n), af(njaf), x(n), soln(nblack),&
   &ctol, rrctol


   integer :: ia(n+1), ja(nja), iaf(nblack+1), jaf(njaf),&
   &idiagf(nblack), RBorder(n), n, nja, njaf, nitmax,&
   &nblack, nred, ierr
   !
   !     local variables
   !
   integer :: i, iter, irpnt, ierror
   double precision :: temp, res0, rvkm, rvk, aconj, omega, resmax
   logical conv
   !
   !     work arrays
   !      avk()      nblack
   !      q()        nblack
   !      aq()       nblack
   !      res()      nblack
   !      v()        nblack
   !      total      5*nblack

   double precision, dimension(:), allocatable :: avk, q, aq, res, v

   allocate( avk(nblack), q(nblack), aq(nblack), res(nblack),&
   &v(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdcnjgd) =="


   !     get [a]{x} and calculate {res} -- note: use {soln} as tmp

   irpnt = nblack+1             ! pointer for REDorder
   if (nred < 1) irpnt = 1


   !     get [a]{x}      note: use soln as tmp

   call xmdrbmtv(a, x, soln, ia, ja, RBorder(1), RBorder(irpnt),&
   &n, nja, nblack, nred)

   res0 = 0.0d0
   do i = 1, nblack
      temp = b( RBorder(i) ) - soln(i)
      res(i) = temp
      res0 = res0 + temp * temp
   enddo


   res0 = dsqrt(res0)

   !     save x into soln and use x as temporary(n) array

   do i = 1, nblack
      soln(i) = x( RBorder(i) )
   enddo

   !
   !     iteration loop
   !

   rvkm = 0.0d0
   do iter = 0, nitmax

      call xmdilusl(v, res, af, iaf, jaf, idiagf, njaf, nblack)
   !
   !      conjugate gradient acceleration
   !
      do i = 1, nblack
         x(RBorder(i)) = v(i)
      enddo
      call xmdrbmtv(a, x, avk, ia, ja, RBorder(1), RBorder(irpnt),&
      &n, nja, nblack, nred)

      rvk = 0.0d+0
      do i = 1, nblack
         rvk = rvk + res(i)*v(i)
      enddo

   !       {q} <- {v} and {aq} <- {avk}

      aconj = 0.0d0
      if (iter > 0) aconj = rvk/rvkm

      do i = 1, nblack
         q(i) = v(i) + aconj * q(i)
         aq(i) = avk(i) + aconj * aq(i)
      enddo

      omega = 0.0d+0
      do i = 1, nblack
         omega = omega + q(i)*aq(i)
      enddo

      omega = rvk/omega
      rvkm = rvk

      resmax = 0.0d+0
      temp = 0.0d+0
      conv = .false.

      do i = 1, nblack
         soln(i) = soln(i) + omega*q(i)
         res(i) = res(i) - omega*aq(i)
         temp = dmax1(dabs(q(i)), dabs(omega*q(i)), temp)
         resmax = resmax + res(i) * res(i)
      enddo

      resmax = dsqrt(resmax)
      conv = temp < ctol
      conv = conv.or.(resmax < rrctol*res0)

      if (conv) go to 890

   enddo

   !     no convergence

   ierr = -1
   !     deallocate ( avk, q, aq, res, v )
   return

890 continue
   nitmax = iter+1

   !     deallocate ( avk, q, aq, res, v )
end subroutine xmdcnjgd





subroutine xmdgtred(a, xx, b, ia, ja, REDorder, nja, n, nred)
   !
   !     update red node solution using black node
   !
   implicit none
   double precision :: a(nja), xx(n), b(n)
   integer :: ia(n+1), ja(nja), REDorder(nred), n, nja, nred
   !
   !     local variables
   !
   integer :: ired, kk, iold

   do ired = 1, nred

      iold = REDorder(ired)

   !     assume diagonal is 1

      xx(iold) = b(iold)

      do kk = ia(iold)+1, ia(iold+1)-1    ! skip diagonal
         xx(iold) = xx(iold) - a(kk) * xx(ja(kk))
      enddo

      xx(iold) = xx(iold) / a( ia(iold) )
   enddo

end subroutine xmdgtred










subroutine xmdilusl(temp, b, af, iaf, jaf, idiagf, njaf,&
&nblack)
   !
   !     input: b, af, iaf, jaf, idiagf, njaf, temp, nblack
   !     output: x
   !
   implicit none
   integer :: njaf, iaf(nblack+1), jaf(njaf), idiagf(nblack), nblack
   double precision :: b(nblack), af(njaf), temp(nblack)
   !
   !     local vars
   !
   integer i, k
   !
   !     forward solve:  Lz=b
   !        (L has unit diagonal)
   !
   do i = 1, nblack
      temp(i) = b( i )
   enddo

   do i = 1, nblack
      do k = iaf(i), idiagf(i)-1
         temp(i) = temp(i) - temp( jaf(k) )*af(k)
      enddo
   enddo
   !
   !        back solve Ux = z
   !           (U does not have unit diag)
   !
   do i = nblack, 1, -1
      do k = idiagf(i)+1, iaf(i+1)-1
         temp(i) = temp(i)-temp( jaf(k) )*af(k)
      enddo
      temp(i) = temp(i)/af( idiagf(i) )
   enddo

end subroutine xmdilusl



subroutine xmdmrgl(i, iaf, jaf, idiagf, list, njaf, n, first,&
&level, levptr, levptr1, nlevptr)
   !
   !      level with drop tolerance
   !          input: levptr   level pointer
   !
   implicit none
   integer :: njaf, n, i, iaf(n+1), jaf(njaf), idiagf(n), list(n),&
   &first, level, levptr(nlevptr), levptr1(n), nlevptr
   !
   !     local
   !
   integer :: irow, next, oldlst, ii, nxtlst, ilevel1, ilevel2
   !
   !     first -> first non-zero in row i of L\U
   !     i = nzero in row k of L\U
   !     list(i) position to next nonzero column
   !
   !     list(last) = n+1
   !
   next = first
10 continue
   if (next < i) then
      oldlst = next         ! current position
      nxtlst = list(next)   ! next column number in link list
      irow = next
   !
   !       scan row "row" of U
   !
      do ii = idiagf(irow)+1, iaf(irow+1)-1
100      continue
         ilevel1 = levptr(ii)+levptr1(irow)+1
         ilevel2 = levptr1(jaf(ii))
         ilevel1 = min0(ilevel1,ilevel2)
         if (ilevel1 <= level) then
            if (jaf(ii) < nxtlst) then
   !
   !     add this index to list()
   !
               list(oldlst) = jaf(ii)
               list(jaf(ii)) = nxtlst
               oldlst = jaf(ii)
               levptr1(jaf(ii)) = ilevel1
            elseif (jaf(ii) == nxtlst) then

               levptr1(jaf(ii)) = ilevel1
               oldlst = jaf(ii)
               nxtlst = list(oldlst)

            elseif (jaf(ii) > nxtlst) then
               oldlst = nxtlst
               nxtlst = list(oldlst)
               goto 100
            endif
         endif
      enddo

      next = list(next)
   !
   !       get next element of l/u
   !
      goto 10
   endif

end subroutine xmdmrgl




subroutine xmdnfctr(a, b, ia, ja, nja, n, ierr)
   !
   !     driver for xmdnfac
   !
   !     input:
   !       a()            coefficient matrix
   !       b()            RHS vector
   !       rwork()        temp real*8 work array
   !       msindx()       index array for integer work array iwork()
   !       ia(), ja()     YSMP pointers for a
   !       nja            size of a() and ja()
   !       njaf           size of af() and jaf()
   !       n              total number of nodes
   !       nblack         number of black nodes
   !
   !     output:
   !       iwork()        integer work array
   !       af()           factored coefficient matrix
   !       ierr           error flag
   !
   use xmdcmn
   use xmdmatrix
   implicit none
   integer :: ia(n+1), ja(nja), n, nja, ierr
   double precision :: a(nja), b(n)

   !
   !     local variables
   !
   integer :: ierror

   if (allocated(af)) deallocate (af)
   allocate( af(njaf), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdnfctr) =="

   call xmdnfac(ia, ja, af, n, nja, njaf, a, b, idiagf, iaf, jaf,&
   &nblack, icolour, RBorder, iblackend)

   if (ierr /= 0) then
      write (miunit,200) ierr
      return
   endif
200 format ('  error in xmdnfctr (xmbnfac)'/'    error flag',i8)

end subroutine xmdnfctr




subroutine xmdnfac(ia, ja, af, n, nja, njaf, a,&
&b, idiagf, iaf, jaf, nblack,&
&icolour, RBorder, iblackend)
   !
   implicit none
   integer :: n, nja, njaf, nblack
   integer :: ia(n+1), ja(nja), idiagf(nblack), iaf(nblack+1),&
   &jaf(njaf), icolour(n), RBorder(n), iblackend(n)
   double precision :: a(nja), af(njaf), b(n)
   !
   !
   !       factors PAP^t  (re-ordered matrix)
   !
   !
   !     local variables
   !
   integer :: iold, id, ii, idd, iii, kk, idkk, iblck, idk, ierror
   double precision :: mult, temp, epsmin

   double precision, dimension(:), allocatable :: row
   integer, dimension(:), allocatable :: list

   allocate( row(nblack), list(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdnfac) =="

   epsmin = 1.0d-300

   row(1:nblack) = 0.0d0
   list(1:nblack) = 0
   !
   !     perform red-black elimination and actural numerical factorization
   !
   !   /                 \ /    \   /    \           /                 \ /    \   /    \
   !   |        :        | |    |   |    | eliminate |        :        | |    |   |    |
   !   |   Dr   :  Arb   | | Xr |   | Br | red nodes |   Dr   :  Arb   | | Xr |   | Br |
   !   |        :        | |    |   |    |  exactly  |        :        | |    |   |    |
   !   |--------:--------| |----| = |----|   =>      |--------:--------| |----| = |----|
   !   |        :        | |    |   |    |           |        :        | |    |   |    |
   !   |  Abr   :  Abb   | | Xb |   | Bb |           |   0    :   R    | | Xb |   | C  |
   !   |        :        | |    |   |    |           |        :        | |    |   |    |
   !   \                 / \    /   \    /           \                 / \    /   \    /
   !
   !      R = Abb - Abr Dr^{-1} Arb    (1)
   !      C = Bb - Abr Dr^{-1} Br      (2)
   !
   blacknodes: do iblck = 1, nblack
   !
   !     scatter packed form of row i
   !
      iold = RBorder(iblck)
      id = iblck
      row(id) = a( ia(iold) )

   !       off-diagonals - black nodes

      do ii = ia(iold)+1, iblackend(iold)
         id = icolour(ja(ii))
         row(id) = row(id) + a(ii) ! Abb  in (1)
      enddo
   !
   !     off-diagonals - red nodes
   !              perform red nodes elimination

      do ii = iblackend(iold)+1, ia(iold+1)-1

         id = ja(ii)  ! red node
         do kk = ia(id)+1, ia(id+1)-1
            idk = ja(kk)
            idkk = icolour(idk)
            row(idkk) = row(idkk) - a(ii) * a(kk) / a( ia(id) )  ! Abr Dr^{-1} Arb  in (1)
         enddo

   !      modify {b}       Bb - Abr Dr^{-1} Br`    in (2)

         b(iold) = b(iold) - a(ii)*b(id) / a( ia(id) )

      enddo
   !
   !     set maker
   !
      do ii = iaf(iblck), iaf(iblck+1)-1
         id = jaf(ii)
         list(id) = 1
      enddo
   !
   !     elimination
   !
      do ii = iaf(iblck), idiagf(iblck)-1

         id = jaf(ii)
         mult = row(id)/af(idiagf(id))
         row(id) = mult

         do iii = idiagf(id)+1, iaf(id+1)-1
            idd = jaf(iii)
            if (list(idd) > 0) then
               row(idd) = row(idd) - mult*af(iii)
            endif
         enddo

      enddo

   !     inverse of diagonal to avoid small pivots

      temp = row(iblck)
      af(idiagf(iblck)) = 1.0d0/(temp + epsmin)
   !
   !     gathering full length array
   !
      do ii = iaf(iblck), iaf(iblck+1)-1
         id = jaf(ii)
         af(ii) = row(id)
         row(id) = 0.0d0
         list(id) = 0
      enddo

   enddo blacknodes

   !     deallocate( row, list )
end subroutine xmdnfac






subroutine xmdorthmn(a, b, x, af, soln,&
&ctol, rrctol,&
&ia, ja, iaf, jaf, idiagf, RBorder,&
&nblack, nred, n, nja, njaf, north, nitmax,&
&ierr)


   implicit none
   double precision :: a(nja), b(n), x(n), af(njaf),&
   &soln(nblack), ctol, rrctol


   integer :: ia(n+1), ja(nja), iaf(nblack+1), jaf(njaf),&
   &idiagf(nblack), RBorder(n),&
   &n, nja, njaf, north, nitmax, nblack, nred, ierr
   !
   !     local variables
   !
   !     work arrays
   !        avk()   nblack
   !        aqaq()  north
   !        q()     (north+1)*nblack
   !        aq()    (north+1)*nblack
   !        res()   nblack
   !
   !       total    2*nblack + 2*(north+1)*nblack + north


   integer :: i, k, iter, j, i2, irpnt, ipnt, ierror
   double precision :: temp, res0, aqr, alpha, omega, resmax, avkaq
   logical :: conv

   double precision, dimension(:), allocatable :: avk, aqaq, q,&
   &aq, res

   allocate( avk(nblack), aqaq(north),&
   &q((north+1)*nblack), aq((north+1)*nblack),&
   &res(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdorthmn) =="

   !
   !     intilaize all vectors
   !
   k = (north+1)*nblack
   do i = 1, k
      q(i) = 0.0d0
      aq(i) = 0.0d0
   enddo

   do i = 1, nblack
   !       v(i) = 0.0d0
      avk(i) = 0.0d0
   enddo

   aqaq(1:north) = 0.0d0

   avk = 0.0d0

   irpnt = nblack+1             ! pointer for REDorder
   if (nred < 1) irpnt = 1

   !     get [a]{x} and calculate {res} -- note: use {soln} as tmp

   call xmdrbmtv(a, x, soln, ia, ja, RBorder(1), RBorder(irpnt),&
   &n, nja, nblack, nred)

   res0 = 0.0d0
   do i = 1, nblack
      temp = b( RBorder(i) ) - soln(i)
      res(i) = temp
      res0 = res0 + temp * temp
   enddo

   res0 = dsqrt(res0)

   !     save x into soln and use x as temporary(n) array

   do i = 1, nblack
      soln(i) = x( RBorder(i) )
   enddo

   !
   !     iteration loop
   !

   do iter = 0, nitmax

      k = mod(iter,north) + 1
      ipnt = (k-1)*nblack

   !     calculate v and substitute into q(ipnt+i)

      call xmdilusl(q(ipnt+1), res, af, iaf, jaf, idiagf, njaf,&
      &nblack)

      alpha = 0.0d0
   !
   !     product of a and v at k
   !
      do i = 1, nblack
         x(RBorder(i)) = q(ipnt+i)
      enddo

      call xmdrbmtv(a, x, avk, ia, ja, RBorder(1), RBorder(irpnt),&
      &n, nja, nblack, nred)


   !     {q} <- {v} and {aq} <- {avk}
   !        note: we already did {q} <- {v}


      do i = 1, nblack
   !         q(ipnt+i) = v(i)
         aq(ipnt+i) = avk(i)
      enddo

      do j = 1, k-1
         i2 = (j-1)*nblack
         avkaq = 0.0d+0
         do i = 1, nblack
            avkaq = avkaq + avk(i)*aq(i2+i)
         enddo
         alpha = -avkaq/aqaq(j)
         do i = 1, nblack
            q(ipnt+i) = q(ipnt+i) + alpha*q(i2+i)
            aq(ipnt+i) = aq(ipnt+i) + alpha*aq(i2+i)
         enddo
      enddo

      aqr = 0.0d+0
      temp = 0.0d+0
      do i = 1, nblack
         aqr = aqr + aq(ipnt+i)*res(i)
         temp = temp + aq(ipnt+i)*aq(ipnt+i)
      enddo
      aqaq(k) = temp

      omega = aqr/aqaq(k)

      resmax = 0.0d+0
      temp = 0.0d+0
      conv = .false.

      do i = 1, nblack
         k = RBorder(i)
         soln(i) = soln(i) + omega*q(ipnt+i)
         res(i) = res(i) - omega*aq(ipnt+i)
         temp = dmax1(dabs(omega*q(ipnt+i)), dabs(q(ipnt+i)), temp)
         resmax = resmax + res(i)*res(i)
      enddo

      resmax = dsqrt(resmax)

      conv = temp < ctol
      conv = conv.or.(resmax < rrctol*res0)

      if (conv) go to 890
   enddo

   !     no convergence

   ierr = -1
   return

890 continue
   nitmax = iter+1
   !     deallocate ( avk, aqaq, q, aq, res )

end subroutine xmdorthmn






subroutine xmdrbmtv(a, xx, amltx, ia, ja, RBorder, REDorder,&
&n, nja, nblack, nred)

   !
   !      multiply [a] and {x} in reduced form
   !
   !      input:
   !             ia(), ja()      YSMP pointers
   !             a()             matrix
   !             xx()            solution vector
   !             n               number of unknowns
   !             nblack          number of black nodes
   !             nred            number of red nodes
   !             nja             size of ja()
   !             RBorder(1:nblack) black node pointer
   !                 original_node_number = RBorder(black_node_number)
   !             REDorder(1:nred)   red node pointer
   !                 original_node_number = REDorder(red_node_number)
   !
   !      output:
   !             amltx()         mulitiplication of [a]{x}in reduced system, i.e., [Abb']{Xb}
   !
   implicit none
   integer :: n, nja, nblack, nred, ia(n+1), ja(nja),&
   &RBorder(nblack), REDorder(nred)
   double precision :: a(nja), amltx(nblack), xx(n)
   !
   !     /                 \ /    \
   !     |        :        | |    |
   !     |   Dr   :  Arb   | | Xr |  after elimination red nodes
   !     |        :        | |    |
   !     |--------:--------| |----|  Abb' = Abb_{diag} + Abb_{off_diag} - Abr * Arb
   !     |        :        | |    |
   !     |  Abr   :  Abb   | | Xb |
   !     |        :        | |    |
   !     \                 / \    /
   !
   !      let's Xr = -Dr^{-1} * Arb * Xb then
   !
   !      Abb' * Xb = Abb_{diag} * Xb + { Abb_{off_diag} - Abr * Dr^{-1} * Arb } Xb
   !                = Abb_{diag} * Xb + ( Abb Abr )_{off_diag} * { Xr+Xb }
   !
   !     local variables
   !
   integer iold, kk, iblck, ired

   !     calculate Xr = -Dr^{-1} * Arb * Xb

   do ired = 1, nred

      iold = REDorder(ired)
      xx(iold) = 0.0d0  ! use red xx(iold) as tmp.

      do kk = ia(iold)+1, ia(iold+1)-1
         xx(iold) = xx(iold) - a(kk) * xx(ja(kk))
      enddo

      xx(iold) = xx(iold) / a( ia(iold) )

   enddo

   !     calculate   Abb' * Xb
   !                  = Abb_{diag} * Xb + ( Abb Abr )_{off_diag} * { Xr+Xb }

   do iblck = 1, nblack

      iold = RBorder(iblck)

   !     obatain Abb_{diag} * Xb
      amltx(iblck) = xx(iold) * a( ia(iold) )

   !     calculate  ( Abb Abr )_{off_diag} * { Xr+Xb } and add to the above
      do kk = ia(iold)+1, ia(iold+1)-1
         amltx(iblck) = amltx(iblck) + a(kk)* xx( ja(kk) )
      enddo

   enddo

end subroutine xmdrbmtv




   !              -------------------------------
   !                   Column Number Rearrange
   !              -------------------------------
   !
   !                      xmdrowrg.f
   !                      version 1.1
   !             last revision: March 29, 2011
   !
   !                (c) 2011   M. Ibaraki
   !
   !     This program rearrange the rows of a matrix so that the diagonal
   !     of each row appear as the first element of each row and
   !     sort row entries in increasing column order.
   !
   !
   !
   !     variable definitions:
   !
   !     ia,ja      YSMP ia,ja data structure
   !     n          total number of unknowns
   !     nja        size of ja array
   !     ierr       error flag
   !
   !     input      ia, n, nja
   !     in-output  ja
   !     output     ierr
   !
   !     required subroutine: xmd_v1.3/xmdshell
   !

subroutine xmdrowrg(ia, ja, n, nja, ierr)

   use xmdcmn
   implicit none
   integer :: n, nja, ierr, ia(n+1), ja(nja)
   !
   !     local variables
   !
   integer :: i, j, node, nentry
   logical :: found

   ierr = 0
   do i = 1, n
      found = .false.
      do j = ia(i), ia(i+1)-1
         if (ja(j) == i) then
            node = ja(ia(i))
            ja(ia(i)) = ja(j)
            ja(j) = node
            found = .true.
         endif
      enddo
      if (.not.found) then
         write (miunit,*) 'error in data structure!!'
         write (miunit,*) 'on the row of ',i
         write (miunit,*) 'the diagonal of this row is missing'
         ierr = 3
         return
      endif
      nentry = ia(i+1)-ia(i)-1
      call xmdshell(ja(ia(i)+1),nentry)
   enddo

end subroutine xmdrowrg


   !      -------------------------------------------------------------
   !        Sparse Matrix Solver for Symmetric and UNSymmetric Matrix
   !      -------------------------------------------------------------
   !
   !                       xmdsolv.f
   !                      version 2.0
   !              last revision:  Feb 22, 2011
   !
   !             (c)copyright 2011   M. Ibaraki
   !
   !     This program solves  [a] x = b  using preconditioned
   !     conjugate gradient type methods.
   !     current version has:
   !                          conjugate gradient
   !                          ORTHOMIN
   !                          CGSTAB
   !
   !     The sparsity pattern of the matrix is stored using ia, ja
   !     data structure.
   !
   !     This subroutine contains acceleration part only
   !
   !
   !     variable definitions:
   !
   !      iacl             choice of acceleration method
   !                       0 = conjugate gradient
   !                       1 = ORTHOMIN
   !                       2 = CGSTAB
   !      a                matrix stored as linear array
   !      af               matrix factors (each row of af
   !                       contains a row L\U)
   !                       where A = LU
   !      aq               A * q^{k}
   !      avk              A * v^{k}
   !      b                right hand side
   !      ctol             convergence tolerance
   !      ia,ja            usual ia, ja arrays
   !      conv             flag for convergence check
   !      iaf, jaf         "ia, ja" arrays for af
   !      idiagf           points to diagonal of af array,
   !                       i.e. af( idiag( i) ) is the
   !                       diagonal element of U, assuming
   !                       A = LU, where L is unit lower triangular
   !      north            number of orthogonalizaion
   !
   !      invord           inverse of lorder:  invord( lorder(i) ) = i
   !
   !      lorder           ordering vector: lorder( new_order ) = old_order
   !      n                number of unknowns
   !      nitmax           max numvber of iterations
   !      nja              size of ja, a, arrays
   !      njaf             size of jaf, af arrays
   !      q                search vector
   !      aqaq             product of aq and aq
   !      res              residual
   !      rwork            temporary n-length vector wkspace
   !      v                (LU)^{-1} r^{k}
   !      x                solution (original ordering)
   !
   !      input        : a,b,ia,ja,iaf,jaf,idiagf,af,n,nja,njaf,lorder
   !                     ctol,rrctol,north
   !      input-output : x,nitmax
   !      temp         : avk,aqaq,q,aq,res,v,rwork
   !
   !
   !
   !     for xmdilusl:
   !       input: res,af,iaf,jaf,idiagf,n,njaf,rwork,lorder
   !              (rwork is a temp vector)
   !       output: v
   !
   !     for xmdcnjgd:
   !       input: ctol,rrctol,a,v,ia,ja,nja,n,iter,invord,lorder
   !       temps: avk
   !       input-output: q,aq,res,x,rvkm
   !       output: conv
   !
   !
   !
   !     for xmdorthmnmn:
   !       input: ctol,rrctol,a,v,ia,ja,nja,n,iter,invord,lorder,north
   !       temps: avk
   !       input-output: q,aq,res,x,aqaq
   !       output: conv
   !

subroutine xmdsolv(a, b, x, ctol, rrctol, ia, ja,&
&nja, n, north, nitmax, iacl, ierr)

   use xmdcmn
   use xmdmatrix

   implicit none
   integer :: ia(n+1), ja(nja), n, nja, north, nitmax, iacl, ierr

   double precision :: a(nja), b(n), x(n), ctol, rrctol
   !
   !     local variables
   !
   integer :: ierror, i, nred

   double precision, dimension(:), allocatable :: soln

   allocate( soln(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdsolv) =="


   nred = n - nblack
   !
   !     conjugate gradient
   !
   if (iacl == 0) then

      call xmdcnjgd(a, b, x, af, soln,&
      &ctol, rrctol,&
      &ia, ja, iaf, jaf, idiagf,&
      &RBorder, nblack,&
      &nred, n, nja, njaf, nitmax, ierr)
   !
   !     ORTHOMIN
   !
   elseif (iacl == 1) then             !  required size

      call xmdorthmn(a, b, x, af, soln,&
      &ctol, rrctol,&
      &ia, ja, iaf, jaf, idiagf,&
      &RBorder,&
      &nblack, nred, n, nja, njaf, north, nitmax,&
      &ierr)
   !
   !     CGSTAB
   !
   elseif (iacl == 2) then         !  required size

      call xmdcgstb(a, b, x, af, soln,&
      &ctol, rrctol,&
      &ia, ja, iaf, jaf, idiagf,&
      &RBorder, nblack,&
      &nred, n, nja, njaf, nitmax, north, ierr)

   endif

   if (ierr == -1) then
   !        write (*,*) 'too many iterations!!'
      ierr = 0
   endif

   !
   !     scatter solution
   !
   do i = 1, nblack               !rgn
      x( RBorder(i) ) = soln(i)    !rgn
   end do                         !rgn
   !       x( RBorder(1:nblack) ) = soln(1:nblack)

   !     recover red nodes
   if (nred > 0) then
      call xmdgtred(a, x, b, ia, ja, RBorder(nblack+1), nja, n, nred)
   endif
   !     deallocate ( soln )
end subroutine xmdsolv



subroutine xmdmrgd(a, af, row, epsrn, i, nblack, ia, iaf, jaf,&
&idiagf, list, RBorder, nja, njaf, n, first,&
&level, levptr, levptr1, nlevptr)

   !
   !      level with drop tolerance
   !
   implicit none
   double precision :: a(nja), af(njaf), row(nblack), epsrn
   integer :: nja, njaf, n, i, nlevptr,&
   &ia(n+1), iaf(nblack+1), jaf(njaf),&
   &idiagf(nblack), list(nblack), RBorder(n),&
   &first, level, levptr(nlevptr), levptr1(nblack), nblack
   !
   !     local
   !
   integer :: irow, next, oldlst, ii, nxtlst, id1, id2,&
   &ilevel1, ilevel2
   double precision mult, tmp1, tmp2
   !
   !     first -> first non-zero in row i of L\U
   !     i = nzero in row k of L\U
   !     list(i) position to next nonzero column
   !
   !     list(last) = n+1
   !
   next = first
10 continue
   if (next < i) then
      oldlst = next         ! current position
      nxtlst = list(next)   ! next column number in link list
      irow = next
   !
   !       scan row "row" of U
   !
      mult = row(irow)/af(idiagf(irow))
      row(irow) = mult

      do ii = idiagf(irow)+1, iaf(irow+1)-1
100      continue
         ilevel1 = levptr(ii)+levptr1(irow)+1
         ilevel2 = levptr1(jaf(ii))
         ilevel1 = min0(ilevel1,ilevel2)
         if (ilevel1 <= level) then
            if (jaf(ii) < nxtlst) then

               id1 = ia(RBorder(i))
               id2 = ia(RBorder(jaf(ii)))
               tmp1 = dabs(mult*af(ii))
               tmp2 = epsrn * dsqrt(dabs(a(id1)*a(id2)))

               if (tmp1 > tmp2) then
   !
   !     add this index to list()
   !
                  list(oldlst) = jaf(ii)
                  list(jaf(ii)) = nxtlst
                  oldlst = jaf(ii)
                  levptr1(jaf(ii)) = ilevel1
                  row(jaf(ii)) = row(jaf(ii)) - mult*af(ii)
               endif
            elseif (jaf(ii) == nxtlst) then

               levptr1(jaf(ii)) = ilevel1
               oldlst = jaf(ii)
               nxtlst = list(oldlst)

               row(jaf(ii)) = row(jaf(ii)) - mult*af(ii)
            elseif (jaf(ii) > nxtlst) then
               oldlst = nxtlst
               nxtlst = list(oldlst)
               goto 100
            endif
         endif
      enddo

      next = list(next)
   !
   !       get next element of l/u
   !
      goto 10
   endif

end subroutine xmdmrgd











subroutine xmdprpc(ia, ja, nja, nn, norder, ierr, redsys)

   !
   !     assign pointers for arrays and ordering vector
   !
   !     input:
   !       ia(), ja()     YSMP pointers for a
   !       iaf(), jaf()   YSMP pointers for af
   !       liwrk          size of integer temp work array iwork()
   !       nja            size of ja()
   !       njaf           size of jaf()
   !       nn             total number of nodes
   !       norder         switch for ordering scheme
   !                        = 0; original ordering
   !                        = 1; RCM ordering
   !                        = 2; Minimum Degree ordering
   !       redsys         flag for reduced system
   !                        = .true.  perform red-black ordering
   !                        = .false. full system treatment
   !     output:
   !       lorder()       ordering vector: lorder( new_order ) = old_order
   !       nblack         number of black nodes
   !       ierr           error flag
   !
   use xmdcmn
   use xmdmatrix

   implicit none
   integer :: ia(nn+1), ja(nja), nja, nn, norder, ierr
   logical :: redsys
   !
   !     local variables
   !
   integer ierror


   allocate(icolour(nn), RBorder(nn), iblackend(nn), lorder(nn),&
   &stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdprpc) =="


   !     row number rearrange

   call xmdrowrg(ia, ja, nn, nja, ierr)

   if (ierr /= 0) then
      write (miunit,150) ierr
      return
   endif
   !
   !     set ordering
   !
   call xmdordng(ia, ja, lorder, nn, nja, norder, ierr)
   if (ierr /= 0) then
      write (miunit,100) ierr
      return
   endif
   !
   !     get red-black ordering
   !

   call xmdRedBlack(ia, ja, lorder, icolour, RBorder,&
   &iblackend, nn, nja, nblack, ierr, redsys)

   allocate(iaf(nblack+1), idiagf(nblack), stat = ierror )
   if (ierror /= 0) stop "== not enough memory (xmdprpc - af) =="

100 format ('  error in xmdprpc (xmdordng)'/'    error flag',i8)
150 format ('  error in xmdprpc (xmdrowrg)'/'    error flag',i8)

end subroutine xmdprpc




subroutine xmdprecl(ia, ja, level, nja, nn, ierr)
   !
   !     assign arrays (ordering, and red-black ordering)
   !     and perform symbolic factorization
   !
   !
   !     input:
   !       ia(), ja()     YSMP pointers for a
   !       iaf(), jaf()   YSMP pointers for af
   !       idiagf()       diagonal pointer for af
   !       iwork()        integer temp work array
   !       level          level of ILU
   !       liwrk          size of integer temp work array iwork()
   !       nblack         number of black nodes
   !       nja            size of ja()
   !       njaf           size of jaf()
   !       njaf           size of af()
   !       nn             total number of nodes
   !       redsys         flag for reduced system
   !                        = .true.  perform red-black ordering
   !                        = .false. full system treatment
   !     output:
   !       ierr           error flag
   !

   use xmdcmn
   use xmdSymFacLev
   use xmdmatrix

   implicit none
   integer :: ia(nn+1), ja(nja), level, nja, nn
   !
   !     local variable
   !
   integer :: ierr
   ierr = 0

   call xmdsfacl(iaf, jaf, idiagf, ia, ja,&
   &icolour, RBorder, iblackend,&
   &nn, nja, njaf, level, nblack, ierr)


   if (ierr /= 0) then
      write (miunit,300) ierr
      return
   endif

300 format ('  error in xmdprecl (xmdsfacl)'/'    error flag',i8)

end subroutine xmdprecl







subroutine xmdprecd(a, b, epsrn, ia, ja, nja, nn, level, ierr)
   !
   !
   !     assign arrays (ordering, and red-black ordering)
   !     and perform symbolic factorization
   !
   !
   !     input:
   !       af()           factored coefficient matrix
   !                      -- note: this array is used as integer work array
   !                               with the size of 2*njaf
   !                      for xmdsfac
   !       ia(), ja()     YSMP pointers for a
   !       iaf(), jaf()   YSMP pointers for af
   !       idiagf()       diagonal pointer for af
   !       iwork()        integer temp work array
   !       level          level of ILU
   !       liwrk          size of integer temp work array iwork()
   !       nblack         number of black nodes
   !       nja            size of ja()
   !       njaf           size of jaf()
   !       njaf           size of af()
   !       nn             total number of nodes
   !     output:
   !       ierr           error flag
   !
   use xmdcmn
   use xmdmatrix
   use xmdSymFacDrop

   implicit none
   double precision :: a(nja), b(nn), epsrn
   integer :: ia(nn+1), ja(nja), level, nja, nn
   !
   !     local variables
   !
   integer :: ierr
   ierr = 0
   !
   !   jaf, af will be calculated
   !
   call xmdsfacd(a, b, af, epsrn, iaf, jaf, idiagf, ia, ja,&
   &icolour, RBorder, iblackend, nn, nja, njaf, level,&
   &nblack, ierr)

   !     print *, 'jaf size ', size(jaf)

   if (ierr /= 0) then
      write (miunit,300) ierr
      return
   endif

300 format ('  error in xmdprecd (xmdsfacd)'/'    error flag',i8)

end subroutine xmdprecd





subroutine xmdcheck(ia, ja, n, nja, ierr)

   !     check the sizes of arrays
   !
   use xmdcmn
   use xmdmatrix

   implicit none
   integer :: ia(n+1), ja(nja), n, nja, ierr
   !
   !     local variables
   !
   integer :: i, iold, id, ii, idd, kk, idkk, iblck, idk, ibk, k

   ierr = 0
   !
   !     check ia, ja data structure
   !
   !
   !     check diag entry in the first column
   !
   do ibk = 1, n
      i = ja( ia(ibk) )
      if (i /= ibk) then
         write (miunit,110) ibk, i
         ierr = 50
         return
      endif
   enddo

110 format ('error in xmdcheck'/&
   &'  diagonal entry is not placed at the first column'/&
   &'  row number:', i5/&
   &'  column number in the first entry:', i5)

   !
   !     check symmetric structure
   !
   do ibk = 1, n
      do k = ia(ibk), ia(ibk+1)-1

         id = ja(k)
         do kk = ia(id), ia(id+1)-1
            if (ja(kk) == ibk) goto 101
         enddo
         write (miunit,105) ibk, id
         ierr = -2

101      continue
      enddo
   enddo

105 format ('warning in xmdcheck'/&
   &'  ia,ja data structure shows non-symmetric structure'/&
   &'  row number:', i5/&
   &'  column number which does not have symmetric part:',&
   &i5)
   !
   !     check red-black ordering connectoins
   !
   blacknodes: do iblck = 1, nblack
      iold = RBorder(iblck)

      offdiagonals: do ii = ia(iold)+1, ia(iold+1)-1
         id = ja(ii)
         idd = icolour(id)
         if (idd < 0) then  !  red node check
            do kk = ia(id)+1, ia(id+1)-1
               idk = ja(kk)
               idkk = icolour(idk)
               if (idkk < 0) then
                  write (miunit,200) id, idk, idkk ; ierr = 3
                  return
               endif
            enddo
         endif
      enddo offdiagonals

   enddo blacknodes

200 format ('error in xmdcheck'/&
   &'  error in red-black ordering'/&
   &'  red node should have black node connections only'/&
   &'  row number:', i5/&
   &'  column number:', i5/&
   &'  column number attribute:',i5)

end subroutine xmdcheck










   !----- subroutine genrcm
   !***************************************************************
   !***************************************************************
   !********   genrcm ..... general Reverse Cuthill Mckee   *******
   !***************************************************************
   !***************************************************************
   !
   !     purpose - genrcm finds the reverse cuthill-mckee
   !        ordering for a general graph. for each connected
   !        component in the graph, genrcm obtains the ordering
   !        by calling the subroutine rcm.
   !
   !     input parameters -
   !        neqns - number of equations
   !        (xadj, adjncy) - array pair containing the adjacency
   !               structure of the graph of the matrix.
   !
   !     output parameter -
   !        perm - vector that contains the rcm ordering.
   !
   !     working parameters -
   !        mask - is used to mark variables that have been
   !               numbered during the ordering process. it is
   !               initialized to 1, and set to zero as each node
   !               is numbered.
   !        xls - the index vector for a level structure.  the
   !               level structure is stored in the currently
   !               unused spaces in the permutation vector perm.
   !
   !     program subroutines -
   !        fnroot, rcm.
   !
   !***************************************************************
   !
subroutine  genrcm ( neqns, nja, xadj, adjncy, perm,&
&mask, xls )
   !
   !***************************************************************
   !
   integer nja, lperm, neqns
   integer adjncy(nja), mask(neqns), perm(neqns),&
   &xls(neqns+1)
   integer xadj(neqns+1), ccsize, i, nlvl, num, root
   !
   !***************************************************************
   !
   do 100 i = 1, neqns
      mask(i) = 1
100 continue
   num = 1
   do 200 i = 1, neqns
   !           ---------------------------------------
   !           for each masked connected component ...
   !           ---------------------------------------
      if (mask(i) .eq. 0) go to 200
      root = i
   !              -----------------------------------------
   !              first find a pseudo-peripheral node root.
   !              note that the level structure found by
   !              fnroot is stored starting at perm(num).
   !              then rcm is called to order the component
   !              using root as the starting node.
   !              -----------------------------------------
   !mi
      lperm = neqns - num +1
   !mi
      call  fnroot ( lperm, neqns, nja, root, xadj, adjncy,&
      &mask, nlvl, xls, perm(num) )
      call     rcm ( lperm, neqns, nja, root, xadj, adjncy,&
      &mask, perm(num), ccsize, xls )
      num = num + ccsize
      if (num .gt. neqns) return
200 continue
   return
end
   !----- subroutine fnroot
   !***************************************************************
   !***************************************************************
   !*******     fnroot ..... find pseudo-peripheral node    *******
   !***************************************************************
   !***************************************************************
   !
   !    purpose - fnroot implements a modified version of the
   !       scheme by gibbs, poole, and stockmeyer to find pseudo-
   !       peripheral nodes.  it determines such a node for the
   !       section subgraph specified by mask and root.
   !
   !    input parameters -
   !       (xadj, adjncy) - adjacency structure pair for the graph.
   !       mask - specifies a section subgraph. nodes for which
   !              mask is zero are ignored by fnroot.
   !
   !    updated parameter -
   !       root - on input, it (along with mask) defines the
   !              component for which a pseudo-peripheral node is
   !              to be found. on output, it is the node obtained.
   !
   !    output parameters -
   !       nlvl - is the number of levels in the level structure
   !              rooted at the node root.
   !       (xls,ls) - the level structure array pair containing
   !                  the level structure found.
   !
   !    program subroutines -
   !       rootls.
   !
   !***************************************************************
   !
subroutine  fnroot ( lls, neqns, nja, root, xadj, adjncy, mask,&
&nlvl, xls, ls )
   !
   !***************************************************************
   !
   integer neqns, nja, lls
   integer adjncy(nja), ls(lls), mask(neqns), xls(neqns+1)
   integer xadj(neqns+1), ccsize, j, jstrt, k, kstop, kstrt,&
   &mindeg, nabor, ndeg, nlvl, node, nunlvl,&
   &root
   !
   !***************************************************************
   !
   !        ---------------------------------------------
   !        determine the level structure rooted at root.
   !        ---------------------------------------------
   call  rootls ( lls, neqns, nja,&
   &root, xadj, adjncy, mask, nlvl, xls, ls )
   ccsize = xls(nlvl+1) - 1
   if ( nlvl .eq. 1 .or. nlvl .eq. ccsize ) return
   !        ----------------------------------------------------
   !        pick a node with minimum degree from the last level.
   !        ----------------------------------------------------
100 jstrt = xls(nlvl)
   mindeg = ccsize
   root = ls(jstrt)
   if ( ccsize .eq. jstrt )  go to 400
   do 300 j = jstrt, ccsize
      node = ls(j)
      ndeg = 0
      kstrt = xadj(node)
      kstop = xadj(node+1) - 1
      do 200 k = kstrt, kstop
         nabor = adjncy(k)
         if ( mask(nabor) .gt. 0 )  ndeg = ndeg + 1
200   continue
      if ( ndeg .ge. mindeg ) go to 300
      root = node
      mindeg = ndeg
300 continue
   !        ----------------------------------------
   !        and generate its rooted level structure.
   !        ----------------------------------------
400 call  rootls ( lls, neqns, nja,&
   &root, xadj, adjncy, mask, nunlvl, xls, ls )
   if (nunlvl .le. nlvl)  return
   nlvl = nunlvl
   if ( nlvl .lt. ccsize )  go to 100
   return
end
   !----- subroutine rcm
   !***************************************************************
   !***************************************************************
   !********     rcm ..... reverse cuthill-mckee ordering   *******
   !***************************************************************
   !***************************************************************
   !
   !     purpose - rcm numbers a connected component specified by
   !        mask and root, using the rcm algorithm.
   !        the numbering is to be started at the node root.
   !
   !     input parameters -
   !        root - is the node that defines the connected
   !               component and it is used as the starting
   !               node for the rcm ordering.
   !        (xadj, adjncy) - adjacency structure pair for
   !               the graph.
   !
   !     updated parameters -
   !        mask - only those nodes with nonzero input mask
   !               values are considered by the routine.  the
   !               nodes numbered by rcm will have their
   !               mask values set to zero.
   !
   !     output parameters -
   !        perm - will contain the rcm ordering.
   !        ccsize - is the size of the connected component
   !               that has been numbered by rcm.
   !
   !     working parameter -
   !        deg - is a temporary vector used to hold the degree
   !               of the nodes in the section graph specified
   !               by mask and root.
   !
   !     program subroutines -
   !        degree.
   !
   !***************************************************************
   !
subroutine  rcm ( llperm, neqns, nja, root, xadj, adjncy,&
&mask, perm, ccsize, deg )
   !
   !***************************************************************
   !
   integer neqns, nja, llperm
   integer adjncy(nja), deg(neqns), mask(neqns),&
   &perm(llperm)
   integer xadj(neqns+1), ccsize, fnbr, i, j, jstop,&
   &jstrt, k, l, lbegin, lnbr, lperm,&
   &lvlend, nbr, node, root
   !
   !***************************************************************
   !
   !        -------------------------------------
   !        find the degrees of the nodes in the
   !        component specified by mask and root.
   !        -------------------------------------
   call  degree ( llperm, neqns, nja,&
   &root, xadj, adjncy, mask, deg,&
   &ccsize, perm )
   mask(root) = 0
   if ( ccsize .le. 1 ) return
   lvlend = 0
   lnbr = 1
   !        --------------------------------------------
   !        lbegin and lvlend point to the beginning and
   !        the end of the current level respectively.
   !        --------------------------------------------
100 lbegin = lvlend + 1
   lvlend = lnbr
   do 600 i = lbegin, lvlend
   !           ----------------------------------
   !           for each node in current level ...
   !           ----------------------------------
      node = perm(i)
      jstrt = xadj(node)
      jstop = xadj(node+1) - 1
   !           ------------------------------------------------
   !           find the unnumbered neighbors of node.
   !           fnbr and lnbr point to the first and last
   !           unnumbered neighbors respectively of the current
   !           node in perm.
   !           ------------------------------------------------
      fnbr = lnbr + 1
      do 200 j = jstrt, jstop
         nbr = adjncy(j)
         if ( mask(nbr) .eq. 0 )  go to 200
         lnbr = lnbr + 1
         mask(nbr) = 0
         perm(lnbr) = nbr
200   continue
      if ( fnbr .ge. lnbr )  go to 600
   !              ------------------------------------------
   !              sort the neighbors of node in increasing
   !              order by degree. linear insertion is used.
   !              ------------------------------------------
      k = fnbr
300   l = k
      k = k + 1
      nbr = perm(k)
400   if ( l .lt. fnbr )  go to 500
      lperm = perm(l)
      if ( deg(lperm) .le. deg(nbr) )  go to 500
      perm(l+1) = lperm
      l = l - 1
      go to 400
500   perm(l+1) = nbr
      if ( k .lt. lnbr )  go to 300
600 continue
   if (lnbr .gt. lvlend) go to 100
   !        ---------------------------------------
   !        we now have the cuthill mckee ordering.
   !        reverse it below ...
   !        ---------------------------------------
   k = ccsize/2
   l = ccsize
   do 700 i = 1, k
      lperm = perm(l)
      perm(l) = perm(i)
      perm(i) = lperm
      l = l - 1
700 continue
   return
end
   !----- subroutine degree
   !***************************************************************
   !***************************************************************
   !********     degree ..... degree in masked component   ********
   !***************************************************************
   !***************************************************************
   !
   !     purpose - this routine computes the degrees of the nodes
   !        in the connected component specified by mask and root.
   !        nodes for which mask is zero are ignored.
   !
   !     input parameter -
   !        root - is the input node that defines the component.
   !        (xadj, adjncy) - adjacency structure pair.
   !        mask - specifies a section subgraph.
   !
   !     output parameters -
   !        deg - array containing the degrees of the nodes in
   !              the component.
   !        ccsize-size of the component specifed by mask and root
   !
   !     working parameter -
   !        ls - a temporary vector used to store the nodes of the
   !               component level by level.
   !
   !***************************************************************
   !
subroutine  degree ( lls, neqns, nja, root, xadj, adjncy, mask,&
&deg, ccsize, ls )
   !
   !***************************************************************
   !
   integer neqns, nja, lls
   integer adjncy(nja), deg(neqns), ls(lls), mask(neqns)
   integer xadj(neqns+1), ccsize, i, ideg, j, jstop, jstrt,&
   &lbegin, lvlend, lvsize, nbr, node, root
   !
   !***************************************************************
   !
   !        -------------------------------------------------
   !        initialization ...
   !        the array xadj is used as a temporary marker to
   !        indicate which nodes have been considered so far.
   !        -------------------------------------------------
   ls(1) = root
   xadj(root) = -xadj(root)
   lvlend = 0
   ccsize = 1
   !        -----------------------------------------------------
   !        lbegin is the pointer to the beginning of the current
   !        level, and lvlend points to the end of this level.
   !        -----------------------------------------------------
100 lbegin = lvlend + 1
   lvlend = ccsize
   !        -----------------------------------------------
   !        find the degrees of nodes in the current level,
   !        and at the same time, generate the next level.
   !        -----------------------------------------------
   do 400 i = lbegin, lvlend
      node = ls(i)
      jstrt = -xadj(node)
      jstop = iabs(xadj(node + 1)) - 1
      ideg = 0
      if ( jstop .lt. jstrt ) go to 300
      do 200 j = jstrt, jstop
         nbr = adjncy(j)
         if ( mask(nbr) .eq. 0 )  go to  200
         ideg = ideg + 1
         if ( xadj(nbr) .lt. 0 ) go to 200
         xadj(nbr) = -xadj(nbr)
         ccsize = ccsize + 1
         ls(ccsize) = nbr
200   continue
300   deg(node) = ideg
400 continue
   !        ------------------------------------------
   !        compute the current level width.
   !        if it is nonzero , generate another level.
   !        ------------------------------------------
   lvsize = ccsize - lvlend
   if ( lvsize .gt. 0 ) go to 100
   !        ------------------------------------------
   !        reset xadj to its correct sign and return.
   !        ------------------------------------------
   do 500 i = 1, ccsize
      node = ls(i)
      xadj(node) = -xadj(node)
500 continue
   return
end
   !----- subroutine rootls
   !***************************************************************
   !***************************************************************
   !********     rootls ..... rooted level structure      *********
   !***************************************************************
   !***************************************************************
   !
   !     purpose - rootls generates the level structure rooted
   !        at the input node called root. only those nodes for
   !        which mask is nonzero will be considered.
   !
   !     input parameters -
   !        root - the node at which the level structure is to
   !               be rooted.
   !        (xadj, adjncy) - adjacency structure pair for the
   !               given graph.
   !        mask - is used to specify a section subgraph. nodes
   !               with mask(i)=0 are ignored.
   !
   !     output parameters -
   !        nlvl - is the number of levels in the level structure.
   !        (xls, ls) - array pair for the rooted level structure.
   !
   !***************************************************************
   !
subroutine  rootls ( lls, neqns, nja, root, xadj, adjncy, mask,&
&nlvl, xls, ls )
   !
   !***************************************************************
   !
   integer neqns, nja, lls
   integer adjncy(nja), ls(lls), mask(neqns), xls(neqns+1)
   integer xadj(neqns+1), i, j, jstop, jstrt, lbegin,&
   &ccsize, lvlend, lvsize, nbr, nlvl,&
   &node, root
   !
   !***************************************************************
   !
   !        ------------------
   !        initialization ...
   !        ------------------
   mask(root) = 0
   ls(1) = root
   nlvl = 0
   lvlend = 0
   ccsize = 1
   !        -----------------------------------------------------
   !        lbegin is the pointer to the beginning of the current
   !        level, and lvlend points to the end of this level.
   !        -----------------------------------------------------
200 lbegin = lvlend + 1
   lvlend = ccsize
   nlvl = nlvl + 1
   xls(nlvl) = lbegin
   !        -------------------------------------------------
   !        generate the next level by finding all the masked
   !        neighbors of nodes in the current level.
   !        -------------------------------------------------
   do 400 i = lbegin, lvlend
      node = ls(i)
      jstrt = xadj(node)
      jstop = xadj(node + 1) - 1
      if ( jstop .lt. jstrt )  go to 400
      do 300 j = jstrt, jstop
         nbr = adjncy(j)
         if (mask(nbr) .eq. 0) go to 300
         ccsize = ccsize + 1
         ls(ccsize) = nbr
         mask(nbr) = 0
300   continue
400 continue
   !        ------------------------------------------
   !        compute the current level width.
   !        if it is nonzero, generate the next level.
   !        ------------------------------------------
   lvsize = ccsize - lvlend
   if (lvsize .gt. 0 ) go to 200
   !        -------------------------------------------------------
   !        reset mask to one for the nodes in the level structure.
   !        -------------------------------------------------------
   xls(nlvl+1) = lvlend + 1
   do 500 i = 1, ccsize
      node = ls(i)
      mask(node) = 1
500 continue
   return
end
















































subroutine odrv(ia, ja, p, ip, isp, n, nja, nsp, flag)



   !mi   subroutine odrv(n, ia,ja,a, p,ip, nsp,isp, path, flag)
   !
   !                                                                3/12/82
   !***********************************************************************
   !  odrv -- driver for sparse matrix reordering routines
   !***********************************************************************
   !
   !  description
   !
   !    odrv finds a minimum degree ordering of the rows and columns
   !    of a matrix m stored in (ia,ja,a) format (see below).  for the
   !    reordered matrix, the work and storage required to perform
   !    gaussian elimination is (usually) significantly less.
   !
   !    note.. odrv and its subordinate routines have been modified to
   !    compute orderings for general matrices, not necessarily having any
   !    symmetry.  the miminum degree ordering is computed for the
   !    structure of the symmetric matrix  m + m-transpose.
   !    modifications to the original odrv module have been made in
   !    the coding in subroutine mdi, and in the initial comments in
   !    subroutines odrv and md.
   !
   !    if only the nonzero entries in the upper triangle of m are being
   !    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
   !    with the diagonal entries placed first in each row.  this is to
   !    ensure that if m(i,j) will be in the upper triangle of m with
   !    respect to the new ordering, then m(i,j) is stored in row i (and
   !    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
   !    strict lower triangle of m, then m(j,i) is stored in row j (and
   !    thus m(i,j) is not stored).
   !
   !
   !  storage of sparse matrices
   !
   !    the nonzero entries of the matrix m are stored row-by-row in the
   !    array a.  to identify the individual nonzero entries in each row,
   !    we need to know in which column each entry lies.  these column
   !    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
   !    ja(k) = j.  to identify the individual rows, we need to know where
   !    each row starts.  these row pointers are stored in the array ia.
   !    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
   !    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
   !    the first location following the last element in the last row.
   !    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
   !    the nonzero entries in the i-th row are stored consecutively in
   !
   !            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
   !
   !    and the corresponding column indices are stored consecutively in
   !
   !            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
   !
   !    since the coefficient matrix is symmetric, only the nonzero entries
   !    in the upper triangle need be stored.  for example, the matrix
   !
   !             ( 1  0  2  3  0 )
   !             ( 0  4  0  0  0 )
   !         m = ( 2  0  5  6  0 )
   !             ( 3  0  6  7  8 )
   !             ( 0  0  0  8  9 )
   !
   !    could be stored as
   !
   !            - 1  2  3  4  5  6  7  8  9 10 11 12 13
   !         ---+--------------------------------------
   !         ia - 1  4  5  8 12 14
   !         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
   !          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
   !
   !    or (symmetrically) as
   !
   !            - 1  2  3  4  5  6  7  8  9
   !         ---+--------------------------
   !         ia - 1  4  5  7  9 10
   !         ja - 1  3  4  2  3  4  4  5  5
   !          a - 1  2  3  4  5  6  7  8  9          .
   !
   !
   !  parameters
   !
   !    n    - order of the matrix
   !
   !    ia   - integer one-dimensional array containing pointers to delimit
   !           rows in ja and a.  dimension = n+1
   !
   !    ja   - integer one-dimensional array containing the column indices
   !           corresponding to the elements of a.  dimension = number of
   !           nonzero entries in (the upper triangle of) m
   !
   !    a    - real one-dimensional array containing the nonzero entries in
   !           (the upper triangle of) m, stored by rows.  dimension =
   !           number of nonzero entries in (the upper triangle of) m
   !
   !    p    - integer one-dimensional array used to return the permutation
   !           of the rows and columns of m corresponding to the minimum
   !           degree ordering.  dimension = n
   !
   !    ip   - integer one-dimensional array used to return the inverse of
   !           the permutation returned in p.  dimension = n
   !
   !    nsp  - declared dimension of the one-dimensional array isp.  nsp
   !           must be at least  3n+4k,  where k is the number of nonzeroes
   !           in the strict upper triangle of m
   !
   !    isp  - integer one-dimensional array used for working storage.
   !           dimension = nsp
   !
   !    path - integer path specification.  values and their meanings are -
   !             1  find minimum degree ordering only
   !             2  find minimum degree ordering and reorder symmetrically
   !                  stored matrix (used when only the nonzero entries in
   !                  the upper triangle of m are being stored)
   !             3  reorder symmetrically stored matrix as specified by
   !                  input permutation (used when an ordering has already
   !                  been determined and only the nonzero entries in the
   !                  upper triangle of m are being stored)
   !             4  same as 2 but put diagonal entries at start of each row
   !             5  same as 3 but put diagonal entries at start of each row
   !
   !    flag - integer error flag.  values and their meanings are -
   !               0    no errors detected
   !              9n+k  insufficient storage in md
   !             10n+1  insufficient storage in odrv
   !             11n+1  illegal path specification
   !
   !
   !  conversion from real to double precision
   !
   !    change the real declarations in odrv and sro to double precision
   !    declarations.
   !
   !-----------------------------------------------------------------------
   !
   integer :: n, nja, nsp
   integer  ia(n+1), ja(nja),  p(n), ip(n),  isp(nsp),  flag,&
   &v, l, head, mmax, next

   !mi
   !mi   original
   !mi
   !mi   integer  ia(*), ja(*),  p(*), ip(*),  isp(*),  path,  flag,
   !mi  *   v, l, head,  tmp, q, n, nsp, mmax, next
   !mi   double precision  a(*)
   !mi   logical  dflag

   integer path
   !mi
   !mi   set for finding ordering only
   !mi
   path = 1
   !mi
   !
   !----initialize error flag and validate path specification
   flag = 0
   if (path.lt.1 .or. 5.lt.path)  go to 111


   !
   !----allocate storage and find minimum degree ordering
   !mi   if ((path-1) * (path-2) * (path-4) .ne. 0)  go to 1


   mmax = (nsp-n)/2
   v    = 1
   l    = v     +  mmax
   head = l     +  mmax
   next = head  +  n
   if (mmax.lt.n)  go to 110
   !
   call md(n, nja, ia,ja, mmax,isp(v),isp(l), isp(head),p,ip,&
   &isp(v), flag)
   if (flag.ne.0)  go to 100
   !mi
   !mic
   !mic----allocate storage and symmetrically reorder matrix
   !mi   1  if ((path-2) * (path-3) * (path-4) * (path-5) .ne. 0)  go to 2
   !mi        tmp = (nsp+1) -      n
   !mi        q   = tmp     - (ia(n+1)-1)
   !mi        if (q.lt.1)  go to 110
   !mic
   !mi        dflag = path.eq.4 .or. path.eq.5
   !mi        call sro
   !mi     *     (n,  ip,  ia, ja, a,  isp(tmp),  isp(q),  dflag)
   !mi
   !
2  return
   !
   ! ** error -- error detected in md
   !             flag = 9 * n + vi from routine mdi.
   !
100 return
   ! ** error -- insufficient storage
110 flag = 10*n + 1
   return
   ! ** error -- illegal path specified
111 flag = 11*n + 1
   return
end



subroutine md(n, nja, ia,ja, mmax, v,l, head,last,next,&
&mark, flag)
   !
   !***********************************************************************
   !  md -- minimum degree algorithm (based on element model)
   !***********************************************************************
   !
   !  description
   !
   !    md finds a minimum degree ordering of the rows and columns of a
   !    general sparse matrix m stored in (ia,ja,a) format.
   !    when the structure of m is nonsymmetric, the ordering is that
   !    obtained for the symmetric matrix  m + m-transpose.
   !
   !
   !  additional parameters
   !
   !    mmax  - declared dimension of the one-dimensional arrays v and l.
   !           mmax must be at least  n+2k,  where k is the number of
   !           nonzeroes in the strict upper triangle of m
   !
   !    v    - integer one-dimensional work array.  dimension = mmax
   !
   !    l    - integer one-dimensional work array.  dimension = mmax
   !
   !    head - integer one-dimensional work array.  dimension = n
   !
   !    last - integer one-dimensional array used to return the permutation
   !           of the rows and columns of m corresponding to the minimum
   !           degree ordering.  dimension = n
   !
   !    next - integer one-dimensional array used to return the inverse of
   !           the permutation returned in last.  dimension = n
   !
   !    mark - integer one-dimensional work array (may be the same as v).
   !           dimension = n
   !
   !    flag - integer error flag.  values and their meanings are -
   !             0      no errors detected
   !             11n+1  insufficient storage in md
   !
   !
   !  definitions of internal parameters
   !
   !    ---------+---------------------------------------------------------
   !    v(s)     - value field of list entry
   !    ---------+---------------------------------------------------------
   !    l(s)     - link field of list entry  (0 =) end of list)
   !    ---------+---------------------------------------------------------
   !    l(vi)    - pointer to element list of uneliminated vertex vi
   !    ---------+---------------------------------------------------------
   !    l(ej)    - pointer to boundary list of active element ej
   !    ---------+---------------------------------------------------------
   !    head(d)  - vj =) vj head of d-list d
   !             -  0 =) no vertex in d-list d
   !
   !
   !             -                  vi uneliminated vertex
   !             -          vi in ek           -       vi not in ek
   !    ---------+-----------------------------+---------------------------
   !    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
   !             -                             -  0 =) vi tail of d-list
   !    ---------+-----------------------------+---------------------------
   !    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
   !             --vk =) compute degree        - vj =) vj last in d-list
   !             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
   !             -  0 =) do not compute degree -
   !    ---------+-----------------------------+---------------------------
   !    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
   !
   !
   !             -                   vi eliminated vertex
   !             -      ei active element      -           otherwise
   !    ---------+-----------------------------+---------------------------
   !    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
   !             -       to be eliminated      -       to be eliminated
   !    ---------+-----------------------------+---------------------------
   !    last(vi) -  m =) size of ei = m        - undefined
   !    ---------+-----------------------------+---------------------------
   !    mark(vi) - -m =) overlap count of ei   - undefined
   !             -       with ek = m           -
   !             - otherwise nonnegative tag   -
   !             -       .lt. mark(vk)         -
   !
   !-----------------------------------------------------------------------
   !
   integer :: n, nja, mmax
   integer  ia(n+1), ja(nja),  v(mmax), l(mmax),  head(n), last(n),&
   &next(n), mark(n),  flag,  tag, dmin, vk,ek, tail, k

   equivalence  (vk,ek)
   !
   !----initialization
   tag = 0
   call  mdi(n, nja, ia,ja, mmax,v,l, head,last,next,&
   &mark,tag, flag)
   if (flag.ne.0)  return
   !
   k = 0
   dmin = 1
   !
   !----while  k .lt. n  do
1  if (k.ge.n)  go to 4
   !
   !------search for vertex of minimum degree
2  if (head(dmin).gt.0)  go to 3
   dmin = dmin + 1
   go to 2
   !
   !------remove vertex vk of minimum degree from degree list
3  vk = head(dmin)
   head(dmin) = next(vk)
   if (head(dmin).gt.0)  last(head(dmin)) = -dmin
   !
   !------number vertex vk, adjust tag, and tag vk
   k = k+1
   next(vk) = -k
   last(ek) = dmin - 1
   tag = tag + last(ek)
   mark(vk) = tag
   !
   !------form element ek from uneliminated neighbors of vk
   call  mdm(n, mmax, vk,tail, v,l, last,next, mark)
   !
   !------purge inactive elements and do mass elimination
   call  mdp(n, mmax, k,ek,tail, v,l, head,last,next, mark)
   !
   !------update degrees of uneliminated vertices in ek
   call  mdu(n, mmax, ek,dmin, v,l, head,last,next, mark)
   !
   go to 1
   !
   !----generate inverse permutation from permutation
4  do 5 k=1,n
      next(k) = -next(k)
5  last(next(k)) = k
   !
   return
end


subroutine mdi(n, nja, ia,ja, mmax,v,l, head,last,next,&
&mark,tag, flag)
   !
   !***********************************************************************
   !  mdi -- initialization
   !***********************************************************************
   integer :: n, nja, mmax
   integer  ia(n+1), ja(nja),  v(mmax), l(mmax),  head(n), last(n),&
   &next(n),&
   &mark(n), tag,  flag,  sfs, vi,dvi, vj, jmin, jmax, j,&
   &lvk, kmax, k, nextvi
   !
   !----initialize degrees, element lists, and degree lists
   do 1 vi=1,n
      mark(vi) = 1
      l(vi) = 0
1  head(vi) = 0
   sfs = n+1
   !
   !----create nonzero structure
   !----for each nonzero entry a(vi,vj)
   do 6 vi=1,n
      jmin = ia(vi)
      jmax = ia(vi+1) - 1
      if (jmin.gt.jmax)  go to 6
      do 5 j=jmin,jmax
         vj = ja(j)
         if (vj-vi) 2, 5, 4
   !
   !------if a(vi,vj) is in strict lower triangle
   !------check for previous occurrence of a(vj,vi)
2        lvk = vi
         kmax = mark(vi) - 1
         if (kmax .eq. 0) go to 4
         do 3 k=1,kmax
            lvk = l(lvk)
            if (v(lvk).eq.vj) go to 5
3        continue
   !----for unentered entries a(vi,vj)
4        if (sfs.ge.mmax)  go to 101
   !
   !------enter vj in element list for vi
         mark(vi) = mark(vi) + 1
         v(sfs) = vj
         l(sfs) = l(vi)
         l(vi) = sfs
         sfs = sfs+1
   !
   !------enter vi in element list for vj
         mark(vj) = mark(vj) + 1
         v(sfs) = vi
         l(sfs) = l(vj)
         l(vj) = sfs
         sfs = sfs+1
5     continue
6  continue
   !
   !----create degree lists and initialize mark vector
   do 7 vi=1,n
      dvi = mark(vi)
      next(vi) = head(dvi)
      head(dvi) = vi
      last(vi) = -dvi
      nextvi = next(vi)
      if (nextvi.gt.0)  last(nextvi) = vi
7  mark(vi) = tag
   !
   return
   !
   ! ** error-  insufficient storage
101 flag = 9*n + vi
   return
end



subroutine mdm(n, mmax, vk,tail, v,l, last,next, mark)
   !
   !***********************************************************************
   !  mdm -- form element from uneliminated neighbors of vk
   !***********************************************************************
   integer :: mmax, n
   integer  vk, tail, v(mmax), l(mmax), last(n), next(n),&
   &mark(n), tag, s,ls,vs,es, b, lb, vb, blp,blpmax
   equivalence  (vs, es)
   !
   !----initialize tag and list of uneliminated neighbors
   tag = mark(vk)
   tail = vk
   !
   !----for each vertex/element vs/es in element list of vk
   ls = l(vk)
1  s = ls
   if (s.eq.0)  go to 5
   ls = l(s)
   vs = v(s)
   if (next(vs).lt.0)  go to 2
   !
   !------if vs is uneliminated vertex, then tag and append to list of
   !------uneliminated neighbors
   mark(vs) = tag
   l(tail) = s
   tail = s
   go to 4
   !
   !------if es is active element, then ...
   !--------for each vertex vb in boundary list of element es
2  lb = l(es)
   blpmax = last(es)
   do 3 blp=1,blpmax
      b = lb
      lb = l(b)
      vb = v(b)
   !
   !----------if vb is untagged vertex, then tag and append to list of
   !----------uneliminated neighbors
      if (mark(vb).ge.tag)  go to 3
      mark(vb) = tag
      l(tail) = b
      tail = b
3  continue
   !
   !--------mark es inactive
   mark(es) = tag
   !
4  go to 1
   !
   !----terminate list of uneliminated neighbors
5  l(tail) = 0
   !
   return
end


subroutine mdp(n, mmax, k,ek,tail, v,l, head,last,next, mark)
   !
   !***********************************************************************
   !  mdp -- purge inactive elements and do mass elimination
   !***********************************************************************
   integer :: n, mmax
   integer  ek, tail,  v(mmax), l(mmax),  head(n), last(n), next(n),&
   &mark(n),  tag, free, li,vi,lvi,evi, s,ls,es, ilp, ilpmax, k, i
   !
   !----initialize tag
   free = 1 !rgn
   tag = mark(ek)
   !
   !----for each vertex vi in ek
   li = ek
   ilpmax = last(ek)
   if (ilpmax.le.0)  go to 12
   do 11 ilp=1,ilpmax
      i = li
      li = l(i)
      vi = v(li)
   !
   !------remove vi from degree list
      if (last(vi).eq.0)  go to 3
      if (last(vi).gt.0)  go to 1
      head(-last(vi)) = next(vi)
      go to 2
1     next(last(vi)) = next(vi)
2     if (next(vi).gt.0)  last(next(vi)) = last(vi)
   !
   !------remove inactive items from element list of vi
3     ls = vi
4     s = ls
      ls = l(s)
      if (ls.eq.0)  go to 6
      es = v(ls)
      if (mark(es).lt.tag)  go to 5
      free = ls
      l(s) = l(ls)
      ls = s
5     go to 4
   !
   !------if vi is interior vertex, then remove from list and eliminate
6     lvi = l(vi)
      if (lvi.ne.0)  go to 7
      l(i) = l(li)
      li = i
   !
      k = k+1
      next(vi) = -k
      last(ek) = last(ek) - 1
      go to 11
   !
   !------else ...
   !--------classify vertex vi
7     if (l(lvi).ne.0)  go to 9
      evi = v(lvi)
      if (next(evi).ge.0)  go to 9
      if (mark(evi).lt.0)  go to 8
   !
   !----------if vi is prototype vertex, then mark as such, initialize
   !----------overlap count for corresponding element, and move vi to end
   !----------of boundary list
      last(vi) = evi
      mark(evi) = -1
      l(tail) = li
      tail = li
      l(i) = l(li)
      li = i
      go to 10
   !
   !----------else if vi is duplicate vertex, then mark as such and adjust
   !----------overlap count for corresponding element
8     last(vi) = 0
      mark(evi) = mark(evi) - 1
      go to 10
   !
   !----------else mark vi to compute degree
9     last(vi) = -ek
   !
   !--------insert ek in element list of vi
10    v(free) = ek
      l(free) = l(vi)
      l(vi) = free
11 continue
   !
   !----terminate boundary list
12 l(tail) = 0
   !
   return
end



subroutine mdu(n, mmax, ek,dmin, v, l, head, last, next, mark)
   !
   !***********************************************************************
   !  mdu -- update degrees of uneliminated vertices in ek
   !***********************************************************************
   integer :: mmax, n
   integer  ek, dmin, v(mmax), l(mmax),  head(n), last(n), next(n),&
   &mark(n),  tag, vi,evi,dvi, s,vs,es, b,vb, ilp,ilpmax,&
   &blp, blpmax, i

   equivalence  (vs, es)
   !
   !----initialize tag
   tag = mark(ek) - last(ek)
   !
   !----for each vertex vi in ek
   i = ek
   ilpmax = last(ek)
   if (ilpmax.le.0)  go to 11
   do 10 ilp=1,ilpmax
      i = l(i)
      vi = v(i)
      if (last(vi))  1, 10, 8
   !
   !------if vi neither prototype nor duplicate vertex, then merge elements
   !------to compute degree
1     tag = tag + 1
      dvi = last(ek)
   !
   !--------for each vertex/element vs/es in element list of vi
      s = l(vi)
2     s = l(s)
      if (s.eq.0)  go to 9
      vs = v(s)
      if (next(vs).lt.0)  go to 3
   !
   !----------if vs is uneliminated vertex, then tag and adjust degree
      mark(vs) = tag
      dvi = dvi + 1
      go to 5
   !
   !----------if es is active element, then expand
   !------------check for outmatched vertex
3     if (mark(es).lt.0)  go to 6
   !
   !------------for each vertex vb in es
      b = es
      blpmax = last(es)
      do 4 blp=1,blpmax
         b = l(b)
         vb = v(b)
   !
   !--------------if vb is untagged, then tag and adjust degree
         if (mark(vb).ge.tag)  go to 4
         mark(vb) = tag
         dvi = dvi + 1
4     continue
   !
5     go to 2
   !
   !------else if vi is outmatched vertex, then adjust overlaps but do not
   !------compute degree
6     last(vi) = 0
      mark(es) = mark(es) - 1
7     s = l(s)
      if (s.eq.0)  go to 10
      es = v(s)
      if (mark(es).lt.0)  mark(es) = mark(es) - 1
      go to 7
   !
   !------else if vi is prototype vertex, then calculate degree by
   !------inclusion/exclusion and reset overlap count
8     evi = last(vi)
      dvi = last(ek) + last(evi) + mark(evi)
      mark(evi) = 0
   !
   !------insert vi in appropriate degree list
9     next(vi) = head(dvi)
      head(dvi) = vi
      last(vi) = -dvi
      if (next(vi).gt.0)  last(next(vi)) = vi
      if (dvi.lt.dmin)  dmin = dvi
   !
10 continue
   !
11 return
end























