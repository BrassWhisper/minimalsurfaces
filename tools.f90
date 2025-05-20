module tools
contains  
  !! Subroutines utilisées pour le calcul matriciel
  ! Affichage à l'écran d'une matrice de taille NxM à coef réels
  subroutine aff_mat_real(A, N, M)
    implicit none
    integer, intent(in) :: N, M
    real(8), dimension(:,:), intent(in) :: A
    integer :: i, j
    
    do i = 1, N
       do j = 1, M
          write(6, '(F10.3,$)') A(i, j)
       end do
       write(6,*)
    end do
  end subroutine aff_mat_real

  ! Affichage à l'écran d'une matrice de taille NxM à coef entiers
  subroutine aff_mat_int(A, N, M)
    implicit none
    integer, intent(in) :: N, M
    integer, dimension(:,:), intent(in) :: A
    integer :: i, j
    
    do i = 1, N
       do j = 1, M
          write(6, '(I4,$)') A(i, j)
       end do
       write(6,*)
    end do
  end subroutine aff_mat_int

  ! Affiche un vecteur de valeurs réelles
  subroutine aff_vect_real(A, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(:), intent(in) :: A
    integer :: i
    
    do i = 1, N
       write(6,*) A(i)
    end do
  end subroutine aff_vect_real

  ! Affiche un vecteur de valeurs entières
  subroutine aff_vect_int(A, N)
    implicit none
    integer, intent(in) :: N
    integer, dimension(:), intent(in) :: A
    integer :: i
    
    do i = 1, N
       write(6,*) A(i)
    end do
  end subroutine aff_vect_int

  ! Réalise le produit AX
  function prod_vect(A, X, N, M)
    implicit none
    real(8), dimension(:), intent(in) :: X
    real(8), dimension(:,:), intent(in) :: A
    integer, intent(in) :: N, M ! Dimension de la matrice A et des vecteurs
    integer :: i, j
    real(8), dimension(N):: prod_vect

    prod_vect = 0D0
    do i = 1, N
       do j = 1, M
          prod_vect(i) = prod_vect(i) + A(i, j) * X(j)
       end do
    end do
  end function prod_vect

  ! Réalise le produit AB
  function prod_mat(A, B, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A, B
    real(8), dimension(N, N) :: prod_mat
    integer :: i, j, k

    prod_mat = 0D0
    do i = 1, N
       do j = 1, N
          do k = 1, N
             prod_mat(i, j) = prod_mat(i, j) + A(i, k) * B(k, j)
          end do
       end do
    end do
  end function prod_mat

  ! Norme 2 d'un vecteur de RN
  real(8) function norme2(X, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N), intent(in) :: X
    integer :: i

    norme2 = 0
    do i = 1, N
       norme2 = norme2 + X(i)**2
    end do
    norme2 = sqrt(norme2)
  end function norme2

  ! Produit tensoriel de 2 vecteurs de R2
  function tensor2(u, v)
    implicit none
    real(8), dimension(2,2) :: tensor2
    real(8), dimension(2) :: u, v

    tensor2(1, 1) = u(1) * v(1)
    tensor2(1, 2) = u(1) * v(2)
    tensor2(2, 1) = u(2) * v(1)
    tensor2(2, 2) = u(2) * v(2)
  end function tensor2

  ! Produit scalaire usuel de 2 vecteurs de R2
  real(8) function scal2(u, v)
    implicit none
    real(8), dimension(2) :: u, v

    scal2 = u(1) * v(1) + u(2) * v(2)
  end function scal2
  
  ! Calcul l'aire d'un triangle avec les coordonnées de ses 3 sommets
  real(8) function area(Z, con, XY)
    implicit none
    real(8), dimension(:,:), intent(in) :: XY
    real(8), dimension(:), intent(in) :: Z
    integer, dimension(3), intent(in) :: con

    !! coordonnées des 3 sommets
    ! x1 = XY(con(1), 1)
    ! x2 = XY(con(2), 1)
    ! x3 = XY(con(3), 1)
    ! y1 = XY(con(1), 2)
    ! y2 = XY(con(2), 2)
    ! y3 = XY(con(3), 2)
    ! z1 = Z(con(1))
    ! z2 = Z(con(2))
    ! z3 = Z(con(3))
    !! coordonnées des 2 vecteurs
    ! ux = x1 - x2
    ! uy = y1 - y2
    ! uz = z1 - z2
    ! vx = x1 - x3
    ! vy = y1 - y3
    ! vz = z1 - z3
    !! coordonnées du produit vectoriel
    ! px = uy * vz - uz * vy
    ! py = uz * vx - ux * vz
    ! pz = ux * vy - uy * vx
    !! Aire
    ! area = sqrt(px**2 + py**2 + pz**2) / 2D0

    area = sqrt(((XY(con(1), 2) - XY(con(2), 2)) * &
         (Z(con(1)) - Z(con(3))) - &
         (Z(con(1)) - Z(con(2))) * &
         (XY(con(1), 2) - XY(con(3), 2)))**2 + &
         ((Z(con(1)) - Z(con(2))) * &
         (XY(con(1), 1) - XY(con(3), 1)) - &
         (XY(con(1), 1) - XY(con(2), 1)) * &
         (Z(con(1)) - Z(con(3))))**2 + &
         ((XY(con(1), 1) - XY(con(2), 1)) * &
         (XY(con(1), 2) - XY(con(3), 2)) - &
         (XY(con(1), 2) - XY(con(2), 2)) * &
         (XY(con(1), 1) - XY(con(3), 1)))**2)&
         / 2D0
  end function area

  ! Calcule l'aire totale de la surface par l'aire des triangles
  real(8) function full_area(C, XY, U, N)
    implicit none
    integer, dimension(:,:), intent(in) :: C
    real(8), dimension(:,:), intent(in) :: XY
    real(8), dimension(:), intent(in) :: U
    integer, intent(in) :: N
    integer :: i

    full_area = 0D0
    do i = 1, 2 * N**2
       full_area = full_area + area(U, C(i, :), XY)
    end do
  end function full_area
end module tools
