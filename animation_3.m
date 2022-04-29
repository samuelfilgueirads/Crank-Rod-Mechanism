%% Atividade 3 - Topicos em Dinamica das Maquinas

clc
clear all
close all

syms R L q A x y q_dot A_dot q_dotdot A_dotdot t

%% Equacoes de loop
f1 = R*cos(q) + L*cos(A) - x; %%componente horizontal 
f2 = R*sin(q) - L*sin(A) - y; %%componente vertical

soln_x = solve(f1==0,x); %%valor de x igualando f1 a zero
soln_y = solve(f2==0,y); %%valor de y igualando f2 a zero

%%Matriz Jacobiana
J = [diff(f1,x) diff(f1,y); diff(f2,x) diff(f2,y)];

%%Matriz b
b = [-diff(f1,q) -diff(f1,A); -diff(f2,q) -diff(f2,A)];

%% Matriz dos coeficientes de velocidade K
K = inv(J)*b;

Kx1 = K(1,1);
Kx2 = K(1,2);
Ky1 = K(2,1);
Ky2 = K(2,2);

%%Matriz velocidade V
V = K*[q_dot; A_dot];

x_dot = V(1); %%velocidade de x
y_dot = V(2); %%velocidade de y

%% Matriz aceleracao A

L1 = diff(K,q); %%matriz coeficiente de aceleracao em relacao a q
L2 = diff(K,A); %%matriz coeficiente de aceleracao em relacao a A

dev_K = L1*q_dot + L2*A_dot; %%derivada de K em relacao a t

Acc = dev_K*[q_dot; A_dot] + K*[q_dotdot; A_dotdot]; %%aceleracao

x_dotdot = Acc(1); %%aceleracao de x
y_dotdot = Acc(2); %%aceleracao de y


n = input('Movimento uniforme? y/n:','s')

switch n
    case 'y'
        
    t = 0:0.0001:0.1; %%incremento do tempo (s)
    AMP = 0.2;
        
    %%carregando variaveis do sistema de 1GDL
    load('q_MU')
    load('q_dot_plot_MU')
    load('q_dotdot_plot_MU')
    load('A_MU')
    load('A_dot_MU')
    load('A_dotdot_MU')
    load('X_MU')
    load('Y_MU')
    load('L_biela')
    load('R_manivela')
    
    L = L_biela;
    R = R_manivela;
    q = q1;
    q_dot = q_dot_plot1;
    A = A1 + (AMP)*sin(q_dot_plot1.*t);
    
    variavel_x = eval(soln_x);
    variavel_y = eval(soln_y);
    
%     figure
%     plot(t,variavel_x,'k')
%     hold on
%     plot(t,X,'r--')
%     hold on
%     plot(t,variavel_y,'b')
%     hold on
%     plot(t,C_plot,'m--')
%     hold off
%     grid on
%     legend('Posicao X (2 GDL)','Posicao X (1 GDL)','Posicao Y (2 GDL)','Posicao Y (1 GDL)')
%     xlabel('Tempo (s)')
%     ylabel('Posicao (mm)')
%     
%     figure
%     plot(variavel_x,variavel_y)
%     
%     figure
%     plot(t,A)
%     hold on
%     plot(t,A1)
%     hold off
    
    
%     %% Grafico 3D
%     [xData, yData, zData] = prepareSurfaceData( variavel_x, variavel_y, q );
% 
%     % Set up fittype and options.
%     ft = 'linearinterp';
% 
%     % Fit model to data.
%     [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );
% 
%     % Plot fit with data.
%     figure( 'Name', 'untitled fit 1' );
%     h = plot( fitresult, [xData, yData], zData );
%     % legend( h, 'untitled fit 1', 'H_A vs. H_X, H_y', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel('Posicao X (mm)')
%     ylabel('Posicao Y (mm)')
%     zlabel('q (rad)')
%     grid on
%     view( -117.0, 12.4 );
% 
%     figure
%     h = plot( fitresult, [xData, yData], zData, 'Style', 'Contour' );
%     legend( h, 'Grau de Amplitude','Location', 'NorthEast' );
%     % Label axes
%     xlabel('Posicao X (mm)')
%     ylabel('Posicao Y (mm)')
%     grid on
% 
%     figure
%     scatter3(variavel_x, variavel_y, q,'filled')
    

    
    otherwise
     
    t = 0:0.0001:3; %%incremento do tempo (s)
    AMP = 0.2;
        
    %%carregando resultados para o sistema de 1GDL
    load('q_MUV')
    load('q_dot_plot_MUV')
    load('q_dotdot_plot_MUV')
    load('A_MUV')
    load('A_dot_MUV')
    load('A_dotdot_MUV')
    load('X_MUV')
    load('Y_MUV')
    load('L_biela')
    load('R_manivela')
    
    L = L_biela;
    R = R_manivela;
    q = q1;
    q_dot = q_dot_plot1;
    A = A1 + (AMP)*sin(q_dot_plot1.*t);
    
    
    variavel_x = eval(soln_x);
    variavel_y = eval(soln_y);
    
%     figure
%     plot(t,variavel_x,'k')
%     hold on
%     plot(t,X,'r--')
%     hold on
%     plot(t,variavel_y,'b')
%     hold on
%     plot(t,C_plot,'m--')
%     hold off
%     grid on
%     legend('Posicao X (2 GDL)','Posicao X (1 GDL)','Posicao Y (2 GDL)','Posicao Y (1 GDL)')
%     xlabel('Tempo (s)')
%     ylabel('Posicao (mm)')
    
    
%     figure
%     plot(variavel_x,variavel_y)
%     
%     figure
%     plot(t,A)
%     hold on
%     plot(t,A1)
%     hold off
%     
    
%     %% Grafico 3D
% figure
%     [xData, yData, zData] = prepareSurfaceData( variavel_x, variavel_y, q );
% 
%     % Set up fittype and options.
%     ft = 'linearinterp';
% 
%     % Fit model to data.
%     [fitresult, gof] = fit( [xData, yData], zData, ft, 'Normalize', 'on' );

    % Plot fit with data.
%     figure( 'Name', 'untitled fit 1' );
%     h = plot( fitresult, [xData, yData], zData );
%     % legend( h, 'untitled fit 1', 'H_A vs. H_X, H_y', 'Location', 'NorthEast' );
%     % Label axes
%     xlabel('Posicao X (mm)')
%     ylabel('Posicao Y (mm)')
%     zlabel('q (rad)')
%     grid on
%     view( -117.0, 12.4 );

%     figure
%     h = plot( fitresult, [xData, yData], zData, 'Style', 'Contour' );
%     legend( h, 'Grau de Amplitude','Location', 'NorthEast' );
%     % Label axes
%     xlabel('Posicao X (mm)')
%     ylabel('Posicao Y (mm)')
%     grid on
% 
%     figure
%     scatter3(variavel_x, variavel_y, q,'filled')
%     xlabel('Posicao X (mm)')
%     ylabel('Posicao Y (mm)')
%     zlabel('Angulo q (rad)')
    
end


%%Ponto de Interesse - Manivela: Coordenadas
Upm = R; %%projecao do ponto de interesse em U
Vpm = 0; %%projecao do ponto de interesse em V

Upm_CM = 0; %%projecao do ponto de interesse em U (mm)
Vpm_CM = -1.016093; %%projecao do ponto de interesse em V (mm)

Xm = Upm*cos(q)-Vpm*sin(q); %%coordenada do ponto de interesse em X
Ym = Upm*sin(q)+Vpm*cos(q); %%coordenada do ponto de interesse em Y

Xm_CM = Upm_CM*cos(q)-Vpm_CM*sin(q); %%coordenada do ponto de interesse em X
Ym_CM = Upm_CM*sin(q)-Vpm_CM*cos(q); %%coordenada do ponto de interesse em Y

%%Ponto de Interesse - Biela: Coordenadas
Upb = 52.17; %%projecao do ponto de interesse em U (mm)
Vpb = 0.000000133; %%projecao do ponto de interesse em V (mm)

Xb = R*cos(q)+Upb*cos(A)+Vpb*sin(A);
Yb = R*sin(q)-Upb*sin(A)+Vpb*sin(A);

%%Ponto de Interesse - Pistao: Coordenadas
Vp = 0;
Up = 0;

Yp = variavel_y+Vp;
Xp = X+Up;

Vp_CM = 0.000332264; %%ponto de interesse em V (mm)
Up_CM = 13.1125; %%ponto de interesse em U (mm)

Yp_CM = variavel_y+Vp_CM;
Xp_CM = X+Up_CM;

Yp_plot = (t.*Yp)./t; %%vetorizando a variavel constante Yp




%%Animacao

%%Vetores posicao

Ypistao = (t.*Yp)./t; %%transformando em vetor Yp

Ypistao_CM = (t.*Yp_CM)./t; %%transformando em vertor Yp_CM

P1 = [0;0]; %%ponto na junta fixa
P2 = [Xm_CM;Ym_CM];%ponto do centro de massa da manivela
P3 = [Xm;Ym];%%ponto na junta de revolucao de ligacao biela-manivela
P4 = [Xb;Yb]; %%ponto do centro de massa da biela
P5 = [Xp;Ypistao];%%ponto da junta de revolucao pistao-biela
P6 = [Xp_CM;Ypistao_CM];%%ponto do centro de massa do pistao


radii = 0.40;

y = figure;

for i=1:10:length(t)
    
   
    ani = subplot(2,2,[1,2,3,4]);
    P1_ponto = viscircles([P1(1) P1(2)],radii,'Color','k');
    P2_ponto = viscircles([P2(1,i) P2(2,i)],radii,'Color','r');
    P3_ponto = viscircles([P3(1,i) P3(2,i)],radii,'Color','k');
%     P4_ponto = viscircles([P4(1,i) P4(2,i)],radii,'Color','r');
    P5_ponto = viscircles([P5(1,i) P5(2,i)],radii,'Color','k');
    P6_ponto = viscircles([P6(1,i) P6(2,i)],radii);
    
    manivela_barra = line([P1(1) P3(1,i)],[P1(2) P3(2,i)],'LineWidth',2.5);
    biela_barra = line([P3(1,i) P5(1,i)],[P3(2,i) P5(2,i)],'LineWidth',2.5);
    
    P4_ponto = viscircles([P4(1,i) P4(2,i)],radii,'Color','r');
    
    y_pontos_pistao = variavel_y(i)+[11 11 -11 -11];
    x_pontos_pistao = Xp(i)+[-28 28 28 -28];

    axis(ani,'equal');
    set(gca,'XLim',[-120 240],'YLim',[-120 120]);
    str = ['Tempo percorrido:' num2str(t(i)) 's'];
    tempo = text(220,-100,str,'FontSize',9);
    
    
    hold on
    
    pistao = fill(x_pontos_pistao,y_pontos_pistao,[.7 .7 .7],'FaceAlpha',0.01);
    
    hold off
    
    view([90 -90])
    ylabel('Eixo Y');
    xlabel('Eixo X')
    
    pause(0.005);
    
%     saveas(y,['Figura' num2str(i)],'png');
    
            delete(P1_ponto)
        delete(P2_ponto)
        delete(P3_ponto)
%         delete(P4_ponto)
        delete(P5_ponto)
%         delete(P6_ponto)
        delete(manivela_barra)
        delete(biela_barra)
        delete(pistao)
        delete(tempo)
        
end


figure
scatter3(Xb, Yb, q)
hold on
scatter3(Xm_CM, Ym_CM, q)
hold off
xlabel('Axis X (mm)')
ylabel('Axis Y (mm)')
zlabel('Angle q (rad)')
legend('Rod Position (P_b)','Crank Position (P_m)')