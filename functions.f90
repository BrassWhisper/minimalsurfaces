!!! Module contenant les fonctions du bord du domaine ainsi que les fonctions
!   phi_i, les fonctions d'aire, de calcul de A_k et b_k...
module functions
  use tools
  implicit none
  real(8), parameter :: min = -1D0, max = 1D0, pi = 3.14159265358
contains
  ! up    ! [min,max]x{min}
  ! down  : [min,max]x{max}
  ! left  : {min}x[min,max]
  ! right : {max}x[min,max]
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Demi cercle sur up et down, 0 sinon !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! min = -1D0, max = 1D0
  ! Calcule la coté haut du domaine carré
  real(8) function up(x, y)
    real(8), intent(in) :: x, y
    up = sqrt(1D0 - x**2)
  end function up
  
  ! Calcule la coté bas du domaine carré
  real(8) function down(x, y)
    real(8), intent(in) :: x, y
    down = sqrt(1D0 - x**2)
  end function down

  ! Calcule la coté gauche du domaine carré
  real(8) function left(x, y)
    real(8), intent(in) :: x, y
    left = 0D0
  end function left

  ! Calcule la coté doit du domaine carré
  real(8) function right(x, y)
    real(8), intent(in) :: x, y
    right = 0D0
  end function right
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Valeur x sur up et down, 0 sur left et 1 sur right !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! min = 0D0, max = 1D0
  ! ! Calcule la coté haut du domaine carré
  ! real(8) function up(x, y)
  !   real(8), intent(in) :: x, y
  !   up = x
  ! end function up
  
  ! ! Calcule la coté bas du domaine carré 
  ! real(8) function down(x, y)
  !   real(8), intent(in) :: x, y
  !   down = x
  ! end function down

  ! ! Calcule la coté gauche du domaine carré
  ! real(8) function left(x, y)
  !   real(8), intent(in) :: x, y
  !   left = 0D0
  ! end function left

  ! ! Calcule la coté doit du domaine carré
  ! real(8) function right(x, y)
  !   real(8), intent(in) :: x, y
  !   right = 1D0
  ! end function right


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Valeur cos(pi * x) sur up et down, cos(pi * y) sur left et right !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! min = -1D0, max = 1D0
  ! ! Calcule la coté haut du domaine carré
  ! real(8) function up(x, y)
  !   real(8), intent(in) :: x, y
  !   up = cos(pi * x)
  ! end function up
  
  ! ! Calcule la coté bas du domaine carré
  ! real(8) function down(x, y)
  !   real(8), intent(in) :: x, y
  !   down = cos(pi * x)
  ! end function down

  ! ! Calcule la coté gauche du domaine carré
  ! real(8) function left(x, y)
  !   real(8), intent(in) :: x, y
  !   left = cos(pi * y)
  ! end function left

  ! ! Calcule la coté doit du domaine carré
  ! real(8) function right(x, y)
  !   real(8), intent(in) :: x, y
  !   right = cos(pi * y)
  ! end function right

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Valeur cos(pi * x) + sin(pi * y) sur le carré !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! Calcule la coté haut du domaine carré
  ! real(8) function up(x, y)
  !   real(8), intent(in) :: x, y
  !   up = cos(pi * x) + sin(pi * y)
  ! end function up
  
  ! ! Calcule la coté bas du domaine carré
  ! real(8) function down(x, y)
  !   real(8), intent(in) :: x, y
  !   down = cos(pi * x) + sin(pi * y)
  ! end function down

  ! ! Calcule la coté gauche du domaine carré
  ! real(8) function left(x, y)
  !   real(8), intent(in) :: x, y
  !   left = cos(pi * x) + sin(pi * y)
  ! end function left

  ! ! Calcule la coté doit du domaine carré
  ! real(8) function right(x, y)
  !   real(8), intent(in) :: x, y
  !   right = cos(pi * x) + sin(pi * y)
  ! end function right

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Valeur x^2 + y^2 sur le carré !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! Calcule la coté haut du domaine carré
  ! real(8) function up(x, y)
  !   real(8), intent(in) :: x, y
  !   up = x**2 + y**2
  ! end function up
  
  ! ! Calcule la coté bas du domaine carré
  ! real(8) function down(x, y)
  !   real(8), intent(in) :: x, y
  !   down = x**2 + y**2
  ! end function down

  ! ! Calcule la coté gauche du domaine carré
  ! real(8) function left(x, y)
  !   real(8), intent(in) :: x, y
  !   left = x**2 + y**2
  ! end function left

  ! ! Calcule la coté doit du domaine carré
  ! real(8) function right(x, y)
  !   real(8), intent(in) :: x, y
  !   right = x**2 + y**2
  ! end function right

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Valeur 0 sur le carré !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! ! Calcule la coté haut du domaine carré
  ! real(8) function up(x, y)
  !   real(8), intent(in) :: x, y
  !   up = 0D0
  ! end function up
  
  ! ! Calcule la coté bas du domaine carré
  ! real(8) function down(x, y)
  !   real(8), intent(in) :: x, y
  !   down = 0D0
  ! end function down

  ! ! Calcule la coté gauche du domaine carré
  ! real(8) function left(x, y)
  !   real(8), intent(in) :: x, y
  !   left = 0D0
  ! end function left

  ! ! Calcule la coté doit du domaine carré
  ! real(8) function right(x, y)
  !   real(8), intent(in) :: x, y
  !   right = 0D0
  ! end function right

  !! Fonctions phi et d_phi de la base de V_h
  ! Calcule la fonction phi_i
  real(8) function phi_i(x, y, i, XY, h)
    real(8), intent(in) :: x, y, h
    integer, intent(in) :: i
    real(8), dimension(:,:), intent(in) :: XY
    
    if (abs(XY(i, 1) - x) > h .or. abs(XY(i, 2) - y) > h) then
       phi_i = 0D0 ! Extérieur
    else
       if ((XY(i, 1) - x) > 0D0) then
          if ((XY(i, 2) - y) > 0D0) then
             if ((XY(i, 1) - x) > (XY(i, 2) - y)) then
                phi_i = (h + XY(i, 1) - x) / h ! Domaine 3
             else
                phi_i = (h + XY(i, 2) - y) / h ! Domaine 1
             end if
          else
             if (((XY(i, 1) - x) + (y - XY(i, 2))) > h) then
                phi_i = 0D0 ! Extérieur
             else
                phi_i = (h + y - XY(i, 2) + XY(i, 1) - x) / h ! Domaine 6
             end if
          end if
       else
          if ((XY(i, 2) - y) < 0D0) then
             if ((x - XY(i, 1)) > (y - XY(i, 2))) then
                phi_i = (h + x - XY(i, 1)) / h ! Domaine 4
             else
                phi_i = (h + y - XY(i, 2)) / h ! Domaine 5
             end if
          else
             if (((x - XY(i, 1)) + (XY(i, 2) - y)) > h) then
                phi_i = 0D0 ! Extérieur
             else
                phi_i = (h + XY(i, 2) - y + x - XY(i, 1)) / h ! Domaine 2
             end if
          end if
       end if
    end if
  end function phi_i

  ! Calcule le gradient de phi_i
  function d_phi_i(x, y, i, XY, h)
    real(8), dimension(2) :: d_phi_i
    real(8), intent(in) :: x, y, h
    integer, intent(in) :: i
    real(8), dimension(:,:), intent(in) :: XY

    if (abs(XY(i, 1) - x) > h .or. abs(XY(i, 2) - y) > h) then
       d_phi_i = (/0D0, 0D0/) ! Extérieur
    else
       if ((XY(i, 1) - x) > 0D0) then
          if (y - (XY(i, 2)) > 0D0) then
             if ((XY(i, 1) - x) > (y - XY(i, 2))) then
                d_phi_i = (/-1D0 / h, 0D0/) ! Domaine 3
             else
                d_phi_i = (/0D0, 1D0 / h/) ! Domaine 1
             end if
          else
             if (((XY(i, 1) - x) + (XY(i, 2) - y)) > h) then
                d_phi_i = (/0D0, 0D0/) ! Extérieur
             else
                d_phi_i = (/-1D0 / h, -1D0 / h/) ! Domaine 6
             end if
          end if
       else
          if ((y - XY(i, 2)) < 0D0) then
             if ((x - XY(i, 1)) > (XY(i, 2) - y)) then
                d_phi_i = (/1D0 / h, 0D0/) ! Domaine 4
             else
                d_phi_i = (/0D0, -1D0 / h/) ! Domaine 5
             end if
          else
             if (((x - XY(i, 1)) + (y - XY(i, 2))) > h) then
                d_phi_i = (/0D0, 0D0/) ! Extérieur
             else
                d_phi_i = (/1D0 / h, 1D0 / h/) ! Domaine 2
             end if
          end if
       end if
    end if
  end function d_phi_i

  ! Trouve dans quel triangle de la grille se situe (x, y)
  ! Renvoie l'indice correspondant à ce triangle dans la matrice C
  integer function find_triangle(x, y, min, h, N)
    implicit none
    real(8), intent(in) :: x, y, min, h
    integer, intent(in) :: N
    real(8) :: line, column

    line = floor((x - min) / h)
    column = floor((y - min) / h)

    if ((x - line * h + y - column * h - 2 * min) > h) then
       find_triangle = 2 * (N * line + column) + 2
    else
       find_triangle = 2 * (N * line + column) + 1
    end if
  end function find_triangle
  
  ! Calcule le gradient du vecteur u
  function d_u(x, y, U, C, XY, min, h, N)
    implicit none
    real(8), dimension(2) :: d_u
    real(8), intent(in) :: x, y, min, h
    real(8), dimension(:), intent(in) :: U
    real(8), dimension(:,:), intent(in) :: XY
    integer, dimension(:,:), intent(in) :: C
    integer, intent(in) :: N
    integer :: i

    i = find_triangle(x, y, min, h, N)

    d_u = U(C(i, 1)) * d_phi_i(x, y, C(i, 1), XY, h) + &
          U(C(i, 2)) * d_phi_i(x, y, C(i, 2), XY, h) + &
          U(C(i, 3)) * d_phi_i(x, y, C(i, 3), XY, h)
  end function d_u

  ! Fonction Ak_c pour former la matrice d'itération
  ! Calcle la matrice 2x2 à l'intérieur de l'intégrale
  function Ak_c(x, y, U, C, XY, min, h, N)
    implicit none
    real(8), dimension(2,2) :: Ak_c
    real(8), intent(in) :: x, y, min, h
    real(8), dimension(:), intent(in) :: U
    integer, dimension(:,:), intent(in) :: C
    real(8), dimension(:,:), intent(in) :: XY
    integer, intent(in) :: N
    real(8) :: temp
    real(8), dimension(2) :: du
    real(8), dimension(2,2) :: tens

    du = d_u(x, y, U, C, XY, min, h, N)
    temp = 1 / sqrt(1 + du(1)**2 + du(2)**2)
    tens = tensor2(du, du)
    
    Ak_c(1, 1) = temp 
    Ak_c(2, 2) = temp
    temp = temp**3

    Ak_c(1, 1) = Ak_c(1, 1) + temp * tens(1, 1) 
    Ak_c(1, 2) = temp * tens(1, 2)
    Ak_c(2, 1) = temp * tens(2, 1)
    Ak_c(2, 2) = Ak_c(2, 2) + temp * tens(2, 2)

  end function Ak_c

  ! Deuxième fonction pour le calcul de Ak_ij (optimisation du temps de calcul)
  real(8) function Ak_ij(Akc, x, y, XY, h, i, j)
    implicit none
    real(8), dimension(2,2), intent(in) :: Akc
    real(8), intent(in) :: x, y, h
    real(8), dimension(:,:), intent(in) :: XY
    integer, intent(in) :: i, j
     
    Ak_ij = scal2(prod_vect(Akc, d_phi_i(x, y, j, XY, h), 2, 2), &
         d_phi_i(x, y, i, XY, h))
  end function Ak_ij
  
  ! Fonction bk_i pour former le second membre
  function bk_i(x, y, U, C, XY, min, h, N, i)
    implicit none
    real(8) :: bk_i
    real(8), intent(in) :: x, y, min, h
    real(8), dimension(:), intent(in) :: U
    integer, dimension(:,:), intent(in) :: C
    real(8), dimension(:,:), intent(in) :: XY
    integer, intent(in) :: N, i
    real(8), dimension(2) :: du, dphi
    
    du = d_u(x, y, U, C, XY, min, h, N)
    dphi = d_phi_i(x, y, i, XY, h)
    bk_i = -(du(1) * dphi(1) + du(2) * dphi(2)) / sqrt(1 + norme2(du, 2)**2)  
  end function bk_i
  
  ! Calcule une fonction test
  ! Utilisée pour tester l'affichage graphique
  real(8) function graph(x, y)
    real(8), intent(in) :: x, y

    !! Goutte d'eau
    ! real(8) :: lambda, A, t, k, ome
    ! lambda = 2D-1
    ! A = 1D-1
    ! k = 1D-1
    ! ome = 1D-1
    ! t = sqrt(x**2 + y**2)
    ! graph = A * cos(t * 20 / (3 * t + 1))

    !! Fonction carrée
    ! graph = x**2 + y**2

    !! Fonction zéro
    graph = 0D0
  end function graph
end module functions
