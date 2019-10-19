clear
clc
 X = imread("x3.bmp");
 figure
 imshow(X);
 title("X input");
 X = double(X);
 
 Y = imread("y5.bmp");
 figure
 imshow(Y);
 title("Y input");
 Y = double(Y);
 
 eps = 1e-12;
 
 function A_Inv = getNextPsevdoInverse(A, delta)
    rows = rows(A);
    columns=columns(A);
    if(rows > columns)
      A_Inv = inv(A'*A + delta^2*eye(columns))*A';
    else
      A_Inv = A'*inv(A*A' + delta^2*eye(rows));
    endif
endfunction
 
function retval = MoorePenrose (A, eps)
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
 retval = A_Inv2;
 endfunction
 
 
 function retval = Greville(A, eps)
    vector = A(1,:)';
    rows = rows(A);
    columns = columns(A);
    
    if(vector.*vector' < eps)
      A_Inv = vector;
    else
      A_Inv = vector/(vector'*vector);
    endif
    
    for i = 2:rows
      vector = A(i,:)';
      Z = eye(columns) - A_Inv * A(1:i-1,:);
      norm = vector'*Z*vector;
      if(norm < eps)
        Z = A_Inv*A_Inv';
        norm = 1+vector'*Z*vector;
      endif
      A_Inv=A_Inv-Z*vector*vector'*A_Inv/norm;
      A_Inv=[A_Inv,Z*vector/norm];
      endfor
    retval = A_Inv;
 endfunction
 
  V = rand(rows(Y), rows(X));
  
  MoorePenrose_X = MoorePenrose(X,eps);
  MoorePenrose_Z = eye(rows(X)) - X*MoorePenrose_X;
  MoorePenrose_A = Y*MoorePenrose_X + V*MoorePenrose_Z';
  MoorePenrose_Y=MoorePenrose_A*X;
  
  figure
  imshow(uint8(MoorePenrose_Y));
  title("Moore-Penrose formula output");
  
  Greville_X = Greville(X, eps);
  Greville_Z = eye(rows(X)) - X*Greville_X;
  Greville_A = Y*Greville_X + V*Greville_Z';
  Greville_Y=Greville_A*X;
  
  figure
  imshow(uint8(Greville_Y));
  title("Greville formula output");