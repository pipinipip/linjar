clc
clear

%% uppgift 3.1 d
disp('uppgift 3.1d')
A = [2.2 -0.4; 0.8 0.4];
disp("Matris A");
disp(A);
% Skapar matrisen A och visar den.

[e1,a1] = eig(A);
% Använder eig funktionen för att ta ut egenvärden och egenvektorer.

eigenvalue_A_1=a1(1,1);
eigenvalue_A_2=a1(2,2);
%Eftersom eig funktionen ger egenvärderna i en matris så tar vi ut värderna
%ur matrisen.

Eigenvalues=['Egenvärderna för A är ',num2str(eigenvalue_A_1), ' och ',num2str(eigenvalue_A_2)];

disp(Eigenvalues);

%Visar egenvärderna


eigenvector_A_1 = [e1(1,1) ; e1(2,1)];
eigenvector_A_2 = [e1(1,2) ; e1(2,2)];
%Eftersom eig funktionen ger egenvektorna i en matris så tar vi ut
%vektorerna ur matrisen.
Eigenvectors=['Egenvektorerna för A är ',mat2str(eigenvector_A_1 , 5), ' och ',mat2str(eigenvector_A_2 , 5)];

disp(Eigenvectors);
%visar egenvektornerna

%% uppgift 3.1 e
disp('uppgift 3.1e')
test1 = A*eigenvector_A_1;

test_vl = ['Vi testar VL ',mat2str(A,5), '*' ,mat2str(eigenvector_A_1, 5), '=' ,mat2str(test1, 5)];

test2 = eigenvalue_A_1*eigenvector_A_1;

test_hl = ['Vi testar HL ',num2str(eigenvalue_A_1, 5), '*', mat2str(eigenvector_A_1, 5), '=', mat2str(test2, 5)];

disp(test_vl)

disp(test_hl)

%vi räknar ut vänsterled och högerled och visar dessa så att de enkelt går
%att jämnföra.
%% uppgift 3.1 f
disp('uppgift 3.1f')
u=eigenvector_A_1+eigenvector_A_2;
% Ger u ett värde

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

%Plottar funktionerna 

%% uppgift 3.2h

disp('uppgift 3.2h')

p1=[0.5 0.25 0.25; 0.25 0.5 0.25; 0.25 0.25 0.5];
%Skapar markovmatrisen p1

p1disp=['P1 är' ,mat2str(p1 , 5)];

disp(p1disp)

%% uppgift 3.2i

disp('uppgift 3.2i')

v=[0.8 ;0.2; 0];
%skapar den initiala tillståndsvektorn v.

p1v=p1*v;

p1vdisp=['P1v är ', mat2str(p1v, 5)];

disp(p1vdisp)
%visar P1v

%% uppgift 3.2j

disp('uppgift 3.2j')

state_v = zeros(3,8);
% Beräknar och lagrar de åtta tillståndvektorerna
for i = 1:8
state_v(:,i) = p1^i*v;
end
plot(transpose(0:8),transpose([v, state_v]));
%plottar värdet för varje steg (0 till 8)

legend("Blåbär","Jordgubb","Vanilj");
%Ger namn till kurverna.

%% uppgift 3.2k

disp('uppgift 3.2k')

p2=[0.6 0.1 0.2; 0.2 0.6 0.1; 0.2 0.3 0.7];
%Skapar markovmatrisen P2
state_v = zeros(3,8);
% Beräknar och lagrar de åtta tillståndvektorerna
for i = 1:8
state_v(:,i) = p2^i*v;
end

plot(transpose(0:8),transpose([v, state_v]));
%plottar värdet för varje steg (0 till 8)

legend("Blåbär","Jordgubb","Vanilj");
%Ger namn till kurverna.
%% uppgift 3.2l

disp('uppgift 3.2l')

p2_till_8=p2^8*v;
%räknar ut p2^8*v

p2_till_8_disp=['övergångsmatrisen vid glasstilfälle 8 är ', mat2str(p2_till_8, 5) ];

disp(p2_till_8_disp)
%visar värdet p2^8*v

%% uppgift 3.3m

disp('uppgift 3.3m')

[e2,a2] = eig(p1);
% Använder eig funktionen för att ta ut egenvärden och egenvektorer.

eigenvalue_p1_1=a2(1,1);
eigenvalue_p1_2=a2(2,2);
eigenvalue_p1_3=a2(3,3);
%Plockar ut värden ut matrisen för att få ut egenvärderna

Eigenvalues_p1=['Egenvärderna för p1 är ',num2str(eigenvalue_p1_1), ' , ',num2str(eigenvalue_p1_2),' och ', num2str(eigenvalue_p1_3)];

disp(Eigenvalues_p1)
%visar egenvärderna.

%% uppgift 3.3o

disp('uppgift 3.3o')

eigenvector_p1_3 = [e2(1,3) ; e2(2,3); e2(3,3)];
%plockar ut Egenvektorn som passar egenvärdet 1. 

Eigenvectror_p1_3_disp=['Eigenvectorn som tillhör egenvärdet 1 är ', mat2str(eigenvector_p1_3, 5)];

disp(Eigenvectror_p1_3_disp)
%visar egenvektorn som passar egenvärdet 1

%% uppgift 3.3p

disp('uppgift 3.3p')

stat_eigenvector_p1_3=(eigenvector_p1_3*(1/(e2(1,3)+e2(2,3)+e2(3,3))));
%skalar om egenvektorn som passar egenvärdet 1 till matrisen p1 så att
%summan blir 1

stat_eigenvector_p1_3_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(stat_eigenvector_p1_3, 5)];

disp(stat_eigenvector_p1_3_disp)
%visar egenvektorn.

%% uppgift 3.3q

disp('uppgift 3.3q')

[e3,a3] = eig(p2);
%använder eig funktionen för att ta ut egenvärden och egenvektorerna till
%p2
eigenvalue_p2_1=a3(1,1);
eigenvalue_p2_2=a3(2,2);
eigenvalue_p2_3=a3(3,3);
%Tar ut egenvärderna ur matrisen som bildas av eig funktionen.

Eigenvalues_p2=['Eigenvalues för p2 är ',num2str(eigenvalue_p2_1), ' , ',num2str(eigenvalue_p2_2),' och ', num2str(eigenvalue_p2_3)];

disp(Eigenvalues_p2)
%visar egenväderna

eigenvector_p2_1 = [e3(1,1) ; e3(2,1); e3(3,1)];
%plockar ut egenvektorn som passar egenvärdet 1 till matrisen.

Eigenvectror_p2_1_disp=['Eigenvectorn som tillhör egenvärdet 1 är ', mat2str(eigenvector_p2_1, 5)];

disp(Eigenvectror_p2_1_disp)
%visar egenvektorn

stat_eigenvector_p2_1=eigenvector_p2_1*(1/(e3(1,1)+e3(2,1)+e3(3,1)));
%skalar om egenvektorn så att summan av dess element blir 1

stat_eigenvector_p2_1_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(stat_eigenvector_p2_1, 5)];

disp(stat_eigenvector_p2_1_disp)
% visar den omskalade vektorn

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
%skapar markovmatrisen B

B4 = B^4;
%efter test så fick vi B^4 var den matris som inte innehöll några nollor.

B4_disp=['B^4 är en matris med bara positiva värden och är ',mat2str(B4, 5)];

disp(B4_disp)
%visar B^4

%% uppgift 3.4v

disp('uppgift 3.4v')

vb= [0; 0; 1; 0; 0; 0; 0; 0; 0];
%Skapar initialvektorn vb

Beata3= B^3;
Beata3v= Beata3*vb;
Beata3_disp=['sannolikheterna att Beata befinner sig på de olika platserna efter 3 förflyttningar är', mat2str(Beata3v, 5) ];
%Räknar ut B^3*vb för att få ut sannolikheten efter 3 steg
disp(Beata3_disp)
%visar B^3*vb

%% uppgift 3.4w

disp('uppgift 3.4w')

[e4,a4] = eig(B);
%använder eig funktionen för att ta ut egenvärdern och egenvektorerna som
%tillhör B

eigenvalue_B_1=a4(1,1);
eigenvalue_B_2=a4(2,2);
eigenvalue_B_3=a4(3,3);
eigenvalue_B_4=a4(4,4);
eigenvalue_B_5=a4(5,5);
eigenvalue_B_6=a4(6,6);
eigenvalue_B_7=a4(7,7);
eigenvalue_B_8=a4(8,8);
eigenvalue_B_9=a4(9,9);
%plockar ut egenvärderna ur matrisen som fås av eig funktionen.

Eigenvalues_B=['Eigenvalues för B är ',num2str(eigenvalue_B_1), ' , ',num2str(eigenvalue_B_2),' , ', num2str(eigenvalue_B_3),num2str(eigenvalue_B_4), ' , ',num2str(eigenvalue_B_5),' , ', num2str(eigenvalue_B_6),num2str(eigenvalue_B_7), ' , ',num2str(eigenvalue_B_8),' och ', num2str(eigenvalue_B_9)];

disp(Eigenvalues_B)
%visar egenvärderna

eigenvector_B_2 = [e4(1,2) ; e4(2,2); e4(3,2);e4(4,2) ; e4(5,2); e4(6,2); e4(7,2) ; e4(8,2); e4(9,2)];
%plockar ut egenvektorn som tillhör egenvärdet 1 i matrisen

eigenvector_B_2_stat=eigenvector_B_2*(1/(e4(1,2) + e4(2,2)+ e4(3,2)+e4(4,2) + e4(5,2)+ e4(6,2)+ e4(7,2) + e4(8,2)+ e4(9,2)));
%skalar om egenvektorn

stat_eigenvector_B_2_disp=['Ifall vi skalar om egenvektorn som tillhör egenvärdet 1 så att alla element i vectorn summerar till 1 så får vi vektorn', mat2str(eigenvector_B_2_stat, 5)];

disp(stat_eigenvector_B_2_disp)
%visar den omskalade egenvektorn
