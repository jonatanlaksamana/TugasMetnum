

% Step 0: isi matrix operator (matrix  A) dan  matrix jawaban (matrix b ) . 

%soal pertama
mtxA = [58.9 0.03 ; -6.10 5.31 ];
mtxB = [59.2 ; 47.0];
 fprintf("jawaban Pertama");
  gaussELWithScalledPivot(mtxA,mtxB)

%soal ke dua
  mtxC = [ 3.3330, 15920, 10.333;
       2.2220,  16.710, 9.6120;
        -1.5611, 5.1792, -1.6855;
       ];
  mtxD = [ 7953; 0.965; 2.714];
  
  
  fprintf("jawaban kedua");
  gaussELWithScalledPivot(mtxC,mtxD)
  
  
  function x = gaussELWithScalledPivot(A,b)
            n = size(A,1); %ambil length matrix operator
            p = (1:n)';	    %set pivot   
          s = max(abs(A'));   %scaled pivot
          for k = 1:(n-1)
            r = abs(A(p(k),k)/s(p(k)));
            kp = k;
            for i = (k+1):n
              t = abs(A(p(i),k)/s(p(i))); 
              if t > r,  r = t;  kp = i;  end
            end
            l = p(kp);  p(kp) = p(k);  p(k) = l;   % interchange p(kp) and p(k) 
            for i = (k+1):n
              A(p(i),k) = A(p(i),k)/A(p(k),k);
              for j = (k+1):n
                A(p(i),j) = A(p(i),j)-A(p(i),k)*A(p(k),j);
              end
            end
          end


        % Step 2: Forward subsitiution to solve  L*y = b, where 
        %         L(i,j) = 0 for j>i; L(i,i) = 1;  L(i,j) = A(p(i),j) for i>j.

          y = zeros(n,1);        % initialize y to be a column vector
          y(1) = b(p(1));
          for i = 2:n
            y(i) = b(p(i));
            for j = 1:(i-1)
              y(i) = y(i)-A(p(i),j)*y(j);
            end
          end
          % Step 3: Back subsitiution to solve  U*x = y
        %         U(i,j) = A(p(i),j) for j>=i; U(i,j) = 0 for i>j.

          x = zeros(n,1);        % initialize x to be a column vector
          x(n) = y(n)/A(p(n),n);
          for i = (n-1):-1:1
            x(i) = y(i);
            for j = (i+1):n
              x(i) = x(i)-A(p(i),j)*x(j);
            end
            x(i) = x(i)/A(p(i),i);
          end
        
  end 

