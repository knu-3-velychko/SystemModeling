clear
clc
 X = imread("x3.bmp");
 X = double(X);
 
 Y = imread("y5.bmp");
 Y = double(Y);
 
 function A_Inv = getNextPsevdoInverse(A, delta)
    rows = rows(A);
    columns=columns(A);
    if(rows > columns)
      A_Inv = inv(A'*A - delta^2*eye(columns))*A';
    else
      A_Inv = A'*inv(A*A' - delta^2*eye(rows));
    endif
endfunction
 
function retval = MoorePenrose (A)
 delta = 100;
 eps = 1e-12;
 diff = 1;
 
 A_Inv1 = getNextPsevdoInverse(A,delta);
 
 while(diff>eps)
    delta = delta/2;
    A_Inv2 = getNextPsevdoInverse(A,delta);
    diff = norm(A_Inv1-A_Inv2);
    A_Inv1=A_Inv2;
 endwhile  
 retval = A_Inv1;
 endfunction
 
  V = rand(rows(Y), rows(X));
  MoorePenrose_X = Moore_Penrose(X);
  MoorePenrose_Z = eye(rows(X)) - X*MoorePenrose_X';
  MoorePenrose_A = Y*MoorePenrose_X + V*MoorePenrose_Z';
  
  imshow(uint8(MoorePenrose_Y));
