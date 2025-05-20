module method
  use tools
  implicit none
  
contains

  ! Calcul de l'inverse d'une matrice diagonale
  function inv_diag(Mat, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: Mat
    real(8), dimension(N, N) :: inv_diag
    integer :: i

    inv_diag = 0D0
    do i = 1, N
       inv_diag(i, i) = 1 / Mat(i, i)
    end do
  end function inv_diag

  ! Calcule la solution de AX = B avec A triangulaire inférieure
  function descente(A, B, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A
    real(8), dimension(N), intent(in) :: B
    real(8), dimension(N) :: descente
    integer :: i, j, k

    descente = 0D0
    do i = 1, N
       do j = 1, i - 1
          descente(i) = descente(i) + A(i, j) * descente(j)
       end do
       descente(i) = (B(i) - descente(i)) / A(i, i)
    end do
  end function descente

  ! Décompose une matrice A en D - E - F avec
  !     D diagonale
  !     E triangulaire inférieure
  !     F triangulaire supérieure
  subroutine decompose(A, D, E, F, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A
    real(8), dimension(N, N), intent(out) :: D, E, F
    integer :: i,j

    D = 0D0
    E = 0D0
    F = 0D0
    do i = 1, N
       do j = 1, N
          if (i == j) then
             D(i, i) = A(i, i)
          else if (i > j) then
             E(i, j) = -A(i, j)
          else if (i < j) then
             F(i, j) = -A(i, j)
          end if
       end do
    end do
  end subroutine decompose
  
  ! Itération d'une méthode
  function iteration_methode(G, X, B, N)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: G
    real(8), dimension(N), intent(in) :: X, B
    real(8), dimension(N) :: iteration_methode

    iteration_methode = prod_vect(G, X, N, N) + B
  end function iteration_methode
  
  ! Méthode de Jacobi pour un système AX = B
  function jacobi(A, X0, B, N, pre)
    implicit none
    ! Arguments de la méthode
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A
    real(8), dimension(N), intent(in) :: X0, B
    real(8), intent(in) :: pre
    ! Résultat de la méthode
    real(8), dimension(N) :: jacobi
    ! Variables
    integer, parameter :: max_it = 1000
    integer :: i
    real(8), dimension(N) :: Xk, Xk1, C
    real(8), dimension(N, N) :: G, D, E, F

    call decompose(A, D, E, F, N)
    G = prod_mat(inv_diag(D, N), E + F, N)
    C = prod_vect(inv_diag(D, N), B, N, N)
    Xk1 = X0
    do i = 1, max_it
       Xk = Xk1
       Xk1 = prod_vect(G, Xk1, N, N) + C
       if (norme2(Xk1 - Xk, N) < pre) then; exit; end if
    end do
    print*, i   
    jacobi = Xk1
  end function jacobi

  ! Méthode de Gauss-Seidel
  function gaussseidel(A, X0, B, N, pre)
    implicit none
    ! Arguments de la méthode
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A
    real(8), dimension(N), intent(in) :: X0, B
    real(8), intent(in) :: pre
    ! Résultat de la méthode
    real(8), dimension(N) :: gaussseidel
    ! Variables
    integer, parameter :: max_it = 100
    integer :: i
    real(8), dimension(N) :: Xk, Xk1, C
    real(8), dimension(N, N) :: G, D, E, F

    call decompose(A, D, E, F, N)
    Xk1 = X0
    do i = 1, max_it
       Xk = Xk1
       Xk1 = descente(D - E, prod_vect(F, Xk1, N, N) + B, N)
       if (norme2(Xk1 - Xk, N) < pre) then; exit; end if
    end do
    print*, i
    gaussseidel = Xk1
  end function gaussseidel

  ! Méthode de relaxation
  function relaxation(A, X0, B, N, pre, omega)
    implicit none
    ! Arguments de la méthode
    integer, intent(in) :: N
    real(8), dimension(N, N), intent(in) :: A
    real(8), dimension(N), intent(in) :: X0, B
    real(8), intent(in) :: pre, omega
    ! Résultat de la méthode
    real(8), dimension(N) :: relaxation
    ! Variables
    integer, parameter :: max_it = 1000
    integer :: i
    real(8), dimension(N) :: Xk, Xk1, C
    real(8), dimension(N, N) :: G, D, E, F

    call decompose(A, D, E, F, N)
    Xk1 = X0
    do i = 1, max_it
       Xk = Xk1
       Xk1 = descente((1D0 / omega) * D - E,&
            prod_vect(((1 - omega) / omega) * D + F, Xk1, N, N) + B,&
            N)
       if (norme2(Xk1 - Xk, N) < pre) then; exit; end if
    end do
    print*, i
    relaxation = Xk1
  end function relaxation
end module method
