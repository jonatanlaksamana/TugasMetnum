% nama : jonathan laksamana purnomo
% GaussEl
clc
clear all
%soal 1
A = [58.9 0.03 ; -6.10 5.31];
B = [59.2 ; 47.0];


%soal-1
A2 = [3.330 15920 10.333 ; 
    2.222 16.710 9.6120 ;
     -1.5611 5.1792 -1.6855];
B2 = [7953; 0.965 ; 2.714];

fprintf("Ini adalah jawaban pertama");
GaussElm(A,B)
fprintf("Ini adalah jawaban kedua");
GaussElm(A2,B2)

function x = GaussElm(A,B)
    n = size(A,1);
    x = zeros(n,1);
    
    for k=1:n-1
        for i=k:n-1
            m = A(i+1,k)/A(k,k);
            for j=k:n
                A(i+1,j)  =A(i+1,j) - m*A(k,j);
                if abs(A(i+1,j)) < 10^(-3)
                    A(i+1,j) = 0 ;
                end
            end
            B(i+1)  =B(i+1) - m*B(k);
        end    
    end
    %backward subs 
    x(n) = B(n,:)/A(n,n); %cari hasil dari element pertama hasil di store dalam variable x 

    for i = n-1:-1:1
       x(i)=(B(i,:) - A(i,i+1:n)  * x(i+1:n,:) )/A(i,i); 
    end
end

