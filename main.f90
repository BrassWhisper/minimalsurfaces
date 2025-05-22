!!! Fichier principal du programme
! Exécute l'itération de la méthode des éléments finis
program elements
  use tools
  use functions
  use method
  implicit none
  !! Variables
  integer, parameter :: N = 30, max_iter = 20
  real(8), parameter :: pre = 1D-2, threshold = 1D-1
  character(len = 10), parameter :: filename='dataplot01'
  real(8), dimension(:,:), allocatable :: Ak, XY
  real(8), dimension(:), allocatable :: U, delta_U, bk
  integer, dimension(:,:), allocatable :: C
  integer :: k
  real(8) :: h
  
  !! correspondance des indices des points du domaine et les indices des matrices C, XY et U
  ! Omega =
  !      min       x        max
  ! min ( 1  . . .  N*(N+1)+1  )
  !     ( .  . . .      .      )
  !  y  ( .  . . .      .      )
  !     ( .  . . .      .      )
  ! max (N+1 . . . (N+1)*(N+1) )

  ! Alloue l'espace mémoire pour les matrices
  allocate(C(2 * N**2, 3), &
       XY((N + 1)**2, 2), &
       U((N + 1)**2), &
       delta_U((N + 1)**2), &
       Ak((N + 1)**2, (N + 1)**2), &
       bk((N + 1)**2))

  ! Initialisation des variables
  Ak = 0D0
  bk = 0D0
  XY = 0D0
  C = 0
  U = 0D0
  delta_U = 0D0
  h = (max - min) / N

  !! Coeur du programme
  call make_connect(C, N)
  call make_XY(XY, N, min, h)
  call init_U(U, N)
  !call aff_mat_int(C, 2 * N**2, 3)
  !call aff_mat_real(XY, (N + 1)**2, 2
  !call aff_vect_real(U, (N + 1)**2)
  
  print*, "Area at iteration :", 0, full_area(C, XY, U, N)

  do k = 1, max_iter
     call make_Ak_temp(Ak, U, C, XY, min, h, N)
     call make_bk_temp(bk, U, C, XY, min, h, N)
     call make_Ak_bk(Ak, bk, XY, N)
     !delta_U = jacobi(Ak, U, bk, (N + 1)**2, pre, 300)
     delta_U = gaussseidel(Ak, U, bk, (N + 1)**2, pre, 300)
     U = U + delta_U
     if (norme2(delta_U, (N + 1)**2) < threshold) then
        print*, "Threshold reached"
        exit
     end if
     print*, "Area at iteration ", k, full_area(C, XY, U, N)
  end do
  
  ! Ecris les valeurs de U dans un fichier
  call plotgraph_U(U, XY, N, filename)
  
contains
  ! Ecrit les points de U dans un fichier pour générer un graphique
  subroutine plotgraph_U(U, XY, N, filename)
    implicit none
    real(8), dimension(:), intent(in) :: U
    real(8), dimension(:,:), intent(in) :: XY
    character(len = 10), intent(in) :: filename
    integer, intent(in) :: N
    integer :: i
    
    open(unit = 10, file = filename)
    do i = 1, (N + 1)**2
       write(10, *)XY(i, 1), XY(i, 2), U(i)
       if (modulo(i, N+1).eq.0) then
          write(10,*)""
       end if
    end do
    close(10)
    print*, "Data written to ", filename
  end subroutine plotgraph_U

  ! Ecrit les points de la fonction graph dans un fichier
  subroutine plotgraph_test(N, XY)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(:,:), intent(in) :: XY 
    integer :: i
    
    open(unit = 10, file = 'dataplot')
    do i = 1, (N + 1)**2
       write(10, *)XY(i, 1), XY(i, 2), graph(XY(i, 1), XY(i, 2))
       if (modulo(i, N+1).eq.0) then
          write(10,*)""
       end if
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

  ! Créé la matrice d'itération Ak temporaire
  ! avant application des conditions de bord
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
    
    area = h**2 / 2D0 ! Aire d'un triangle
    Ak = 0D0
    ! Boucle d'itération pour tous les triangles
    do t = 1, 2 * N**2, 2
       !! Triangle formé des points C(t,:)
       ! Point à l'intérieur du triangle
       xc = XY(C(t, 1), 1) + (h / 3D0)
       yc = XY(C(t, 1), 2) + (h / 3D0)
       Akc = Ak_c(xc, yc, U, C, XY, min, h, N)
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
  ! avant application des conditions de bord
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
    
    area = h**2 / 2D0 ! Aire d'un triangle
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
  ! Applique les conditions de bord
  subroutine make_Ak_bk(Ak, bk, XY, N)
    real(8), dimension(:,:), intent(inout) :: Ak
    real(8), dimension(:), intent(inout) :: bk
    real(8), dimension(:,:), intent(in) :: XY
    integer, intent(in) :: N
    integer :: i, j
    do i = 1, (N + 1)**2
       if (i < (N + 2) .or. i > ((N + 1)**2 - (N + 1)) .or. &
            modulo(i, N+1).eq.1 .or. modulo(i, N+1).eq.0) then
          do j = 1, (N + 1)**2
             Ak(i, j) = 0D0
          end do
          Ak(i, i) = 1D0
          bk(i) = 0D0
       end if
    end do
  end subroutine make_Ak_bk

  ! Initialise le vecteur U
  ! Valeur aux bord selon les fonctions choisie, 0 sinon
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
    end do
  end subroutine init_U
end program elements
