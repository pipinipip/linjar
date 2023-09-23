clc
clear

%% uppgift 3.1 d
disp('uppgift 3.1d')
A = [2.2 -0.4; 0.8 0.4];
disp("Matris A");
disp(A);

[e1,a1] = eig(A);

eigenvalue_A_1=a1(1,1);
eigenvalue_A_2=a1(2,2);

Eigenvalues=['Eigenvalues för A är ',num2str(eigenvalue_A_1), ' och ',num2str(eigenvalue_A_2)];

disp(Eigenvalues);

eigenvector_A_1 = [e1(1,1) ; e1(2,1)];
eigenvector_A_2 = [e1(1,2) ; e1(2,2)];

Eigenvectors=['Eigenvectors för A är ',mat2str(eigenvector_A_1 , 5), ' och ',mat2str(eigenvector_A_2 , 5)];

disp(Eigenvectors);

%% uppgift 3.1 e
disp('uppgift 3.1e')
test1 = A*eigenvector_A_1;

test_vl = ['Vi testar VL ',mat2str(A,5), '*' ,mat2str(eigenvector_A_1, 5), '=' ,mat2str(test1, 5)];

test2 = eigenvalue_A_1*eigenvector_A_1;

test_hl = ['Vi testar HL ',num2str(eigenvalue_A_1, 5), '*', mat2str(eigenvector_A_1, 5), '=', mat2str(test2, 5)];

disp(test_vl)

disp(test_hl)

%% uppgift 3.1 f
disp('uppgift 3.1f')
u=eigenvector_A_1+eigenvector_A_2;

hold on
plotv(eigenvector_A_1)
plotv(eigenvector_A_2)
plotv(u)
plotv(A*u)
plotv(A^2*u)
xlabel ('x-axeln')
ylabel ('y-axeln')
legend('eigenvector_1=x1','eigenvector_2=x2','u=eigenvector_1+eigenvector_2','Au','A^2u')
hold off

%% uppgift 3.2h

disp('uppgift 3.2h')

p1=[0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5];

p1disp=['P1 är' ,mat2str(p1 , 5)];

disp(p1disp)

%% uppgift 3.2i

disp('uppgift 3.2i')

v=[0.8 ;0.2; 0];

p1v=p1*v;

p1vdisp=['P1v är ', mat2str(p1v, 5)];

disp(p1vdisp)

%% uppgift 3.2j

disp('uppgift 3.2j')

state_v = zeros(3,8);
% Beräkna och lagra de åtta tillståndvektorerna
for i = 1:8
state_v(:,i) = p1^i*v;
end
% Plotta resultatet, med initiala värdena först,
plot(transpose(0:8),transpose([v, state_v]));
% Ange förklarande text i figuren
legend("Blåbär","Jordgubb","Vanilj");

%% uppgift 3.2k

disp('uppgift 3.2k')

p2=[0.6 0.1 0.2; 0.2 0.6 0.1; 0.2 0.3 0.7];

state_v = zeros(3,8);
% Beräkna och lagra de åtta tillståndvektorerna
for i = 1:8
state_v(:,i) = p2^i*v;
end
% Plotta resultatet, med initiala värdena först,
plot(transpose(0:8),transpose([v, state_v]));
% Ange förklarande text i figuren
legend("Blåbär","Jordgubb","Vanilj");

%% uppgift 3.2l

disp('uppgift 3.2l')

p2_till_8=p2^8*v;

p2_till_8_disp=['övergångsmatrisen vid glasstilfälle 8 är ', mat2str(p2_till_8, 5) ];

disp(p2_till_8_disp)

%% uppgift 3.3m

disp('uppgift 3.3m')

[e2,a2] = eig(p1);

eigenvalue_p1_1=a2(1,1);
eigenvalue_p1_2=a2(2,2);
eigenvalue_p1_3=a2(3,3);

Eigenvalues_p1=['Eigenvalues för p1 är ',num2str(eigenvalue_p1_1), ' , ',num2str(eigenvalue_p1_2),' och ', num2str(eigenvalue_p1_3)];

disp(Eigenvalues_p1)

%% uppgift 3.3o

disp('uppgift 3.3o')

eigenvector_p1_3 = [e2(1,3) ; e2(2,3); e2(3,3)];

Eigenvectror_p1_3_disp=['Eigenvectorn som tillhör egenvärdet 1 är ', mat2str(eigenvector_p1_3, 5)];

disp(Eigenvectror_p1_3_disp)

%% uppgift 3.3p

disp('uppgift 3.3p')

stat_eigenvector_p1_3=(eigenvector_p1_3/0.57735)/3;

%fixa detta ^^^^^^ så det är som på uppgift 3.3q

stat_eigenvector_p1_3_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(stat_eigenvector_p1_3, 5)];

disp(stat_eigenvector_p1_3_disp)

%% uppgift 3.3q

disp('uppgift 3.3q')

[e3,a3] = eig(p2);

eigenvalue_p2_1=a3(1,1);
eigenvalue_p2_2=a3(2,2);
eigenvalue_p2_3=a3(3,3);

Eigenvalues_p2=['Eigenvalues för p2 är ',num2str(eigenvalue_p2_1), ' , ',num2str(eigenvalue_p2_2),' och ', num2str(eigenvalue_p2_3)];

disp(Eigenvalues_p2)

eigenvector_p2_1 = [e3(1,1) ; e3(2,1); e3(3,1)];

Eigenvectror_p2_1_disp=['Eigenvectorn som tillhör egenvärdet 1 är ', mat2str(eigenvector_p2_1, 5)];

disp(Eigenvectror_p2_1_disp)

stat_eigenvector_p2_1=eigenvector_p2_1*(1/(e3(1,1)+e3(2,1)+e3(3,1)));

stat_eigenvector_p2_1_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(stat_eigenvector_p2_1, 5)];

disp(stat_eigenvector_p2_1_disp)

%% uppgift 3.4t

disp('uppgift 3.4t')

B=[1/3 1/4 0 1/4 0 0 0 0 0;
 1/3 1/4 1/3 0 1/5 0 0 0 0;
   0 1/4 1/3 0 0 1/4 0 0 0;
 1/3 0 0 1/4 1/5 0 1/3 0 0;
 0 1/4 0 1/4 1/5 1/4 0 1/4 0;
 0 0 1/3 0 1/5 1/4 0 0 1/3;
 0 0 0 1/4 0 0 1/3 1/4 0;
 0 0 0 0 1/5 0 1/3 1/4 1/3;
 0 0 0 0 0 1/4 0 1/4 1/3];

B4 = B^4;

B4_disp=['B^4 är en matris med bara positiva värden och är ',mat2str(B4, 5)];

disp(B4_disp)

%% uppgift 3.4v

disp('uppgift 3.4v')

vb= [0; 0; 1; 0; 0; 0; 0; 0; 0];

Beata3= B^3;
Beata3v= Beata3*vb;
Beata3_disp=['sannolikheterna att Beata befinner sig på de olika platserna efter 3 förflyttningar är', mat2str(Beata3v, 5) ];

disp(Beata3_disp)


%% uppgift 3.4w

disp('uppgift 3.4w')

[e4,a4] = eig(B);

eigenvalue_B_1=a4(1,1);
eigenvalue_B_2=a4(2,2);
eigenvalue_B_3=a4(3,3);
eigenvalue_B_4=a4(4,4);
eigenvalue_B_5=a4(5,5);
eigenvalue_B_6=a4(6,6);
eigenvalue_B_7=a4(7,7);
eigenvalue_B_8=a4(8,8);
eigenvalue_B_9=a4(9,9);

Eigenvalues_B=['Eigenvalues för B är ',num2str(eigenvalue_B_1), ' , ',num2str(eigenvalue_B_2),' , ', num2str(eigenvalue_B_3),num2str(eigenvalue_B_4), ' , ',num2str(eigenvalue_B_5),' , ', num2str(eigenvalue_B_6),num2str(eigenvalue_B_7), ' , ',num2str(eigenvalue_B_8),' och ', num2str(eigenvalue_B_9)];

disp(Eigenvalues_B)

eigenvector_B_2 = [e4(1,2) ; e4(2,2); e4(3,2);e4(4,2) ; e4(5,2); e4(6,2); e4(7,2) ; e4(8,2); e4(9,2)];

eigenvector_B_2_stat=eigenvector_B_2*(1/(e4(1,2) + e4(2,2)+ e4(3,2)+e4(4,2) + e4(5,2)+ e4(6,2)+ e4(7,2) + e4(8,2)+ e4(9,2)));

stat_eigenvector_B_2_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(eigenvector_B_2_stat, 5)];

disp(stat_eigenvector_B_2_disp)
