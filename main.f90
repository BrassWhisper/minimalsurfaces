program elements
  use tools
  use functions
  use method
  implicit none
  !! Variables
  integer, parameter :: N = 50, max_iter = 10
  real(8), parameter :: pre = 1D-2
  real(8), dimension(:,:), allocatable :: Ak, XY
  real(8), dimension(:), allocatable :: Z, U, delta_U, bk
  integer, dimension(:,:), allocatable :: C
  integer :: t, t1, t2, t3, i, k
  real(8) :: x, y, h
  real(8), dimension(2) :: d1, d2, d3
  
  ! Alloue l'espace mémoire pour les matrices
  allocate(C(2 * N**2, 3), &
       XY((N + 1)**2, 2), &
       Z((N + 1)**2), &
       U((N + 1)**2), &
       delta_U((N + 1)**2), &
       Ak((N + 1)**2, (N + 1)**2), &
       bk((N + 1)**2))

  ! Initialisation des variables
  Ak = 0D0
  bk = 0D0
  XY = 0D0
  Z = 0D0
  C = 0
  U = 0D0
  delta_U = 0D0
  h = (max - min) / N

  !! correspondance des indices entre la matrice Ak et les matrices C, XY et Z
  ! Ak =
  ! ( 1  . . .  N*(N+1)+1  )
  ! ( .  . . .      .      )
  ! ( .  . . .      .      )
  ! ( .  . . .      .      )
  ! (N+1 . . . (N+1)*(N+1) )
  ! Ak(i,j)=(i-1)*(N+1)+j

  !! Coeur du programme
  !call plotgraph(N, min, max, Ak)
  call make_connect(C, N)
  call make_XY(XY, N, min, h)
  call make_Z(Z, N, Ak)
  !call aff_mat_int(C, 2 * N**2, 3)
  !call aff_mat_real(XY, (N + 1)**2, 2)
  !call aff_vect_real(Z, (N + 1)**2)
  
  ! do i = 1, (N + 1)**2
  !    U(i) = exp(-5D1 * ((XY(i, 1) - 5D-1)**2 + ((XY(i, 2) - 5D-1)**2)))
  ! end do
  ! print*, full_area(C, XY, U, N)
  
  ! call make_Ak_temp(Ak, U, C, XY, min, h, N)
  ! call aff_mat_real(Ak, (N + 1)**2, (N + 1)**2)
  ! call make_bk_temp(bk, U, C, XY, min, h, N)
  ! call aff_vect_real(bk, (N + 1)**2)
  ! call make_Ak_bk(Ak, bk, XY, N)
  ! !call aff_mat_real(Ak, (N + 1)**2, (N + 1)**2)
  ! call aff_vect_real(bk, (N + 1)**2)

  ! call aff_vect_real(U, (N + 1)**2)
  ! U1 = jacobi(Ak, U, bk, (N + 1)**2, pre)
  ! call aff_vect_real(U + U1, (N + 1)**2)

  ! call aff_vect_real(prod_vect(Ak, U + U1, (N + 1)**2, (N + 1)**2), (N + 1)**2)

  call init_U(U, N)
  print*, full_area(C, XY, U, N)
  call aff_vect_real(U, (N + 1)**2)
  do k = 1, max_iter
     call make_Ak_temp(Ak, U, C, XY, min, h, N)
     call make_bk_temp(bk, U, C, XY, min, h, N)
     call make_Ak_bk(Ak, bk, XY, N)
     !delta_U = jacobi(Ak, U, bk, (N + 1)**2, pre)
     delta_U = gaussseidel(Ak, U, bk, (N + 1)**2, pre)
     U = U + delta_U
     ! print*, full_area(C, XY, U, N)
     ! print*,""
     ! call aff_vect_real(prod_vect(Ak, delta_U, (N + 1)**2, (N + 1)**2), (N + 1)**2)
     ! print*,""
     ! call aff_vect_real(bk, (N + 1)**2)
  end do
  
  print*,""
  call aff_vect_real(U, (N + 1)**2)

  call plotgraph_U(U, XY, N, 'dataplotU3')
  
contains
  ! Ecrit les points de U dans un fichier pour générer un graphique
  subroutine plotgraph_U(U, XY, N, filename)
    implicit none
    real(8), dimension(:), intent(in) :: U
    real(8), dimension(:,:), intent(in) :: XY
    character(len = 10), intent(in) :: filename
    integer, intent(in) :: N
    integer :: i
    real(8) :: x, y, h
    
    open(unit = 10, file = filename)
    do i = 1, (N + 1)**2
       write(10, *)XY(i, 1), XY(i, 2), U(i)
    end do
    close(10)
  end subroutine plotgraph_U

  ! Ecrit les points de la fonction graph dans un fichier
  subroutine plotgraph_test(N, min, max, mat)
    implicit none
    integer, intent(in) :: N
    real(8), intent(in) :: min, max
    real(8), dimension(:,:), intent(inout) :: mat
    integer :: i, j
    real(8) :: x, y, h

    h = (max - min) / N
  
    open(unit = 10, file = 'dataplot')
    do i = 1, N + 1
       x = min + (i - 1) * h
       do j = 1, N + 1
          y = min + (j - 1) * h
          mat(i, j) = graph(x, y)
          write(10, *)x, y, mat(i, j)
       end do
    end do
    close(10)
  end subroutine plotgraph_test

  ! Créé la matrice de connectivité
  subroutine make_connect(C, N)
    implicit none
    integer, intent(in) :: N
    integer, dimension(:,:) :: C
    integer :: i, j, ind_C, ind_mat

    do i = 1, N
       do j = 1, N
          ind_C = (i-1) * N + j
          ind_mat = (i-1) * (N+1) + j
          C(2 * ind_C - 1, 1) = ind_mat
          C(2 * ind_C - 1, 2) = ind_mat + 1
          C(2 * ind_C - 1, 3) = ind_mat + N + 1
          C(2 * ind_C,     1) = ind_mat + N + 2
          C(2 * ind_C,     2) = ind_mat + N + 1
          C(2 * ind_C,     3) = ind_mat + 1
       end do
    end do
  end subroutine make_connect

  ! Créé la matrice des coordonnées x et y
  subroutine make_XY(XY, N, min, h)
    implicit none
    real(8), dimension(:,:), intent(inout) :: XY
    integer, intent(in) :: N
    real(8), intent(in) :: min, h
    integer :: i, j

    do i = 0, N
       do j = 0, N
          XY(i * (N+1) + j + 1, :) = (/h * i + min, h * j + min/)
       end do
    end do
  end subroutine make_XY

  ! Créé la matrice des coordonnées z
  subroutine make_Z(Z, N, mat)
    implicit none
    real(8), dimension(:), intent(inout) :: Z
    real(8), dimension(:,:), intent(in) :: mat
    integer, intent(in) :: N
    integer :: i, j

    do i = 1, N+1
       do j = 1, N+1
          Z((i-1) * (N+1) + j) = mat(i, j)
       end do
    end do
  end subroutine make_Z

  ! Créé la matrice d'itération Ak temporaire
  ! avant application des conditions de bords
  subroutine make_Ak_temp(Ak, U, C, XY, min, h, N)
    implicit none
    real(8), dimension(:,:), intent(out) :: Ak
    real(8), dimension(:), intent(in) :: U
    integer, dimension(:,:), intent(in) :: C
    real(8), dimension(:,:), intent(in) :: XY
    real(8), intent(in) :: min, h
    integer, intent(in) :: N
    integer :: t, t1, i
    real(8) :: area, xc, yc
    real(8), dimension(2,2) :: Akc
    
    area = 1D0 / h**2 ! Aire d'un triangle
    Ak = 0D0
    ! Boucle d'itération pour tous les triangles
    do t = 1, 2 * N**2, 2
       !! Triangle formé des points C(t,:)
       ! Point à l'intérieur du triangle
       xc = XY(C(t, 1), 1) + (h / 3D0)
       yc = XY(C(t, 1), 2) + (h / 3D0)
       Akc = Ak_c(xc, yc, U, C, XY, min, h, N)
       !print*, C(t,:)
       !print*, Akc
       Ak(C(t, 1), C(t, 1)) = Ak(C(t, 1), C(t, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 1), C(t, 1))
       Ak(C(t, 1), C(t, 2)) = Ak(C(t, 1), C(t, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 1), C(t, 2))
       Ak(C(t, 1), C(t, 3)) = Ak(C(t, 1), C(t, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 1), C(t, 3))
       Ak(C(t, 2), C(t, 1)) = Ak(C(t, 2), C(t, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 2), C(t, 1))
       Ak(C(t, 2), C(t, 2)) = Ak(C(t, 2), C(t, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 2), C(t, 2))
       Ak(C(t, 2), C(t, 3)) = Ak(C(t, 2), C(t, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 2), C(t, 3))
       Ak(C(t, 3), C(t, 1)) = Ak(C(t, 3), C(t, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 3), C(t, 1))
       Ak(C(t, 3), C(t, 2)) = Ak(C(t, 3), C(t, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 3), C(t, 2))
       Ak(C(t, 3), C(t, 3)) = Ak(C(t, 3), C(t, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t, 3), C(t, 3))
   
       !! Triangle formé des points C(t+1,:)
       t1 = t + 1
       ! Point à l'intérieur du triangle
       xc = XY(C(t1, 1), 1) - (h / 3D0)
       yc = XY(C(t1, 1), 2) - (h / 3D0)
       Akc = Ak_c(xc, yc, U, C, XY, min, h, N)
       !print*, "C", C(t+1,:)
       !print*, "xc", xc, "yc", yc
       !print*, "Akc", Akc
       !print*, C(t1, 1), Ak_ij(Akc, xc, yc, XY, h, C(t1, 1), C(t1, 1))
       !print*, d_phi_i(xc, yc, C(t1, 1), XY, h)
       !print*, "find", C(find_triangle(xc, yc, min, h, N),:)
       !print*, C(t1, 2), Ak_ij(Akc, xc, yc, XY, h, C(t1, 2), C(t1, 2))
       !print*, C(t1, 3), Ak_ij(Akc, xc, yc, XY, h, C(t1, 3), C(t1, 3))
       Ak(C(t1, 1), C(t1, 1)) = Ak(C(t1, 1), C(t1, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 1), C(t1, 1))
       Ak(C(t1, 1), C(t1, 2)) = Ak(C(t1, 1), C(t1, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 1), C(t1, 2))
       Ak(C(t1, 1), C(t1, 3)) = Ak(C(t1, 1), C(t1, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 1), C(t1, 3))
       Ak(C(t1, 2), C(t1, 1)) = Ak(C(t1, 2), C(t1, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 2), C(t1, 1))
       Ak(C(t1, 2), C(t1, 2)) = Ak(C(t1, 2), C(t1, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 2), C(t1, 2))
       Ak(C(t1, 2), C(t1, 3)) = Ak(C(t1, 2), C(t1, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 2), C(t1, 3))
       Ak(C(t1, 3), C(t1, 1)) = Ak(C(t1, 3), C(t1, 1)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 3), C(t1, 1))
       Ak(C(t1, 3), C(t1, 2)) = Ak(C(t1, 3), C(t1, 2)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 3), C(t1, 2))
       Ak(C(t1, 3), C(t1, 3)) = Ak(C(t1, 3), C(t1, 3)) + &
            area * Ak_ij(Akc, xc, yc, XY, h, C(t1, 3), C(t1, 3))
       end do
  end subroutine make_Ak_temp

  ! Créé le vecteur bk temporaire
  ! avant application des condition de bord
  subroutine make_bk_temp(bk, U, C, XY, min, h, N)
    implicit none
    real(8), dimension(:), intent(out) :: bk
    real(8), dimension(:), intent(in) :: U
    integer, dimension(:,:), intent(in) :: C
    real(8), dimension(:,:), intent(in) :: XY
    real(8), intent(in) :: min, h
    integer, intent(in) :: N
    integer :: t, t1
    real(8) :: area, xc, yc
    
    area = 1D0 / h**2 ! Aire d'un triangle
    bk = 0D0
    ! Boucle d'itération pour tous les triangles
    do t = 1, 2 * N**2, 2
       !! Triangle formé des points C(t,:)
       ! Point à l'intérieur du triangle
       xc = XY(C(t, 1), 1) + (h / 3D0)
       yc = XY(C(t, 1), 2) + (h / 3D0)
       bk(C(t, 1)) = bk(C(t, 1)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t, 1))
       bk(C(t, 2)) = bk(C(t, 2)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t, 2))
       bk(C(t, 3)) = bk(C(t, 3)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t, 3))

       !! Triangle formé des points C(t+1,:)
       t1 = t + 1
       ! Point à l'intérieur du triangle
       xc = XY(C(t1, 1), 1) - (h / 3D0)
       yc = XY(C(t1, 1), 2) - (h / 3D0)
       bk(C(t1, 1)) = bk(C(t1, 1)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t1, 1))
       bk(C(t1, 2)) = bk(C(t1, 2)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t1, 2))
       bk(C(t1, 3)) = bk(C(t1, 3)) + &
            area * bk_i(xc, yc, U, C, XY, min, h, N, C(t1, 3))
    end do
  end subroutine make_bk_temp

  ! Créé la matrice Ak et le vecteur bk
  subroutine make_Ak_bk(Ak, bk, XY, N)
    real(8), dimension(:,:), intent(inout) :: Ak
    real(8), dimension(:), intent(inout) :: bk
    real(8), dimension(:,:), intent(in) :: XY !,Ak_temp
    !real(8), dimension(:), intent(in) :: bk_temp
    integer, intent(in) :: N
    integer :: i, j!, k, f, m
    !real(8), dimension(4*N,4*N) :: Ak_fixed
    !real(8), dimension(4*N) :: bk_fixed

    ! k = 1
    ! f = 1
    do i = 1, (N + 1)**2
       if (i < (N + 2) .or. i > ((N + 1)**2 - (N + 1)) .or. &
            modulo(i, N+1).eq.1 .or. modulo(i, N+1).eq.0) then
          !! Dans le cas d'une réduction de la taille de Ak
          ! ! Génère la matrice des conditions au bord
          ! m = 0
          ! do j = f, 4 * N
          !    Ak_fixed(f, j) = Ak_temp(i, i + m)
          !    Ak_fixed(j, f) = Ak_temp(i + m, i)
          !    m = m + 1
          ! end do
          ! ! Génère le vecteur de conditions au bord
          ! bk_fixed(f) = bk_temp(i)
          ! f = f + 1
          !! Dans le cas ou on garde Ak de la même taille
          do j = 1, (N + 1)**2
             Ak(i, j) = 0D0
          end do
          Ak(i, i) = 1D0
          ! if (i < (N + 2)) then
          !    bk(i) = left(XY(i, 1), XY(i, 2))
          ! else
          !    if (i > ((N + 1)**2 - (N + 1))) then
          !       bk(i) = right(XY(i, 1), XY(i, 2))
          !    else
          !       if (modulo(i, N+1).eq.1) then
          !          bk(i) = up(XY(i, 1), XY(i, 2))
          !       else
          !          if (modulo(i, N+1).eq.0) then
          !             bk(i) = down(XY(i, 1), XY(i, 2))
          !          end if
          !       end if
          !    end if
          ! end if
          bk(i) = 0D0 
       !else
          !! Dans le cas d'une réduction de la taille de Ak
          ! ! Génère la matrice d'itération
          ! m = 0
          ! do j = k, (N + 1)**2 - 4 * N
          !    Ak(k, j) = Ak_temp(i, i + m)
          !    Ak(j, k) = Ak_temp(i + m, i)
          !    m = m + 1
          ! end do
          ! ! Génère le vecteur des itérations
          ! bk(k) = bk_temp(i)
          ! k = k + 1
       end if
    end do
  end subroutine make_Ak_bk

  ! Initialise le vecteur U
  subroutine init_U(U, N)
    implicit none
    real(8), dimension(:), intent(out) :: U
    integer, intent(in) :: N
    integer :: i

    do i = 1, (N + 1)**2
       if (i < (N + 2)) then
          U(i) = left(XY(i, 1), XY(i, 2))
       else
          if (i > ((N + 1)**2 - (N + 1))) then
             U(i) = right(XY(i, 1), XY(i, 2))
          else
             if (modulo(i, N+1).eq.1) then
                U(i) = up(XY(i, 1), XY(i, 2))
             else
                if (modulo(i, N+1).eq.0) then
                   U(i) = down(XY(i, 1), XY(i, 2))
                else
                   U(i) = 0D0
                end if
             end if
          end if
       end if
       print*, i, U(i)
  end do
  end subroutine init_U
end program elements
