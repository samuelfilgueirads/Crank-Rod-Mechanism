clc
close all
clear all

%%Declaracao dos parametros importantes
L = 1170.2; %%mm
R = 400.1; %%mm
C = 0; %%mm

t = 0:0.0001:0.024; %%incremento do tempo (s)

n = input('Uniform motion? y/n:','s')

switch n
    case 'y'
        
q_dotdot = 0; %%aceleracao angular da generalizada (rad/s^2)
q_dotdot_plot = (t.*q_dotdot)./t; %%transformar aceleracao angular da generalizada como vetor de mesmo numero de elementos de t para plotagem mais adiante
q_dot = 5000*((2*pi)/60); %%velocidade angular da generalizada (rad/s)
q_dot_plot = (t.*q_dot)./t; %%transformar velocidade angular da generalizada como vetor de mesmo numero de elementos de t para plotagem mais adiante
q_zero = 0; %%angulo inicial da generalizada (rad)
q = q_zero + t.*q_dot; %%angulo da generalizada (rad) %%movimento uniforme: s = s0+vt
 
    otherwise

q_dotdot = 20; %%aceleracao angular da generalizada (rad/s^2)
q_dotdot_plot = (t.*q_dotdot)./t; %%transformar aceleracao angular da generalizada como vetor de mesmo numero de elementos de t para plotagem mais adiante
q_dot_zero = 0; %%velocidade angular inicial da generalizada (rad/s)
q_dot = q_dot_zero + t.*q_dotdot; %%velocidade angular da generalizada (rad/s) %%movimento uniformemente variado v = v0 + at
q_dot_plot = q_dot; %%para plotagem da velocidade angular da generalizada
q_zero = 0; %%angulo inicial da generalizada (rad)
q = q_zero + t.*q_dot_zero + (q_dotdot*(t.^2)./2); %%angulo da generalizada (rad) %%movimento uniformemente variado: s = s0+v0t+(a(t^2)/2)
     
end

%%Analise Cinematica

A = asin((R*sin(q) - C)/L); %%angulo A (rad)
X = R*cos(q) + L*cos(A); %%posicao x (definir medida)

Ka = R*cos(q)./(L*cos(A)); %%coeficiente de velocidade de A
Kx = (-R*L*sin(q).*cos(A) - R*L*cos(q).*sin(A))./(L*cos(A)); %%coeficiente de velocidade de X

A_dot = Ka.*q_dot; %%velocidade de A
X_dot = Kx.*q_dot; %%velocidade de X

La = ((-R*L*sin(q).*cos(A)) + (R*L*cos(q).*sin(A).*Ka))./((L*cos(A)).^2); %%coeficiente de aceleracao de A
Lx = -R*cos(q) - L*La.*sin(A) - L*(Ka.^2).*cos(A); %%coeficiente de aceleracao de X

A_dotdot = Ka.*q_dotdot + La.*(q_dot.^2); %%aceleracao de A
X_dotdot = Kx.*q_dotdot + Lx.*(q_dot.^2); %%aceleracao de X

%%Posicao: Pontos de Interesse

%%Ponto de Interesse - Manivela: Coordenadas
Upm = R; %%projecao do ponto de interesse em U
Vpm = 0; %%projecao do ponto de interesse em V

Xm = Upm*cos(q)-Vpm*sin(q); %%coordenada do ponto de interesse em X
Ym = Upm*sin(q)-Vpm*cos(q); %%coordenada do ponto de interesse em Y

%%Ponto de Interesse - Biela: Coordenadas
Upb = L/2; %%projecao do ponto de interesse em U
Vpb = 0; %%projecao do ponto de interesse em V

Xb = R*cos(q)+Upb*cos(A)+Vpb*sin(A);
Yb = R*sin(q)-Upb*sin(A)+Vpb*sin(A);

%%Ponto de Interesse - Pistao: Coordenadas
Vp = 0;
Up = 0;

Yp = C+Vp;
Xp = X+Up;

Yp_plot = (t.*Yp)./t; %%vetorizando a variavel constante Yp

%%Velocidade: Pontos de Interesse

%%Manivela
Kxm = (-Upm*sin(q)-Vpm*cos(q)); %%coeficiente de velocidade de Xm
Xm_dot = Kxm.*(q_dot); %%velocidade de Xm

Kym = Upm*cos(q)-Vpm*sin(q);%%coeficiente de velocidade de Ym
Ym_dot = Kym.*(q_dot);%%velocidade de Ym

%%Biela
Kxb = (-R*sin(q)+(-Upb*sin(A)+Vpb*cos(A)).*Ka);%%coeficiente de velocidade de Xb
Xb_dot = Kxb.*(q_dot);%%velocidade de Xb

Kyb = (R*cos(q)+(-Upb*cos(A)-Vpb*sin(A)).*Ka);%%coeficiente de velocidade de Yb
Yb_dot = Kyb.*(q_dot);%%velocidade de Yb

%%Pistao
Yp_dot = 0; %%velocidade de Yp
Xp_dot = X_dot; %%velocidade de Xp


%%Aceleracao: Pontos de Interesse

%%Manivela
Lxm = -Upm*cos(q)+Vpm*sin(q); %%coeficiente de aceleracao de Xm
Xm_dotdot = Kxm.*(q_dotdot)+Lxm.*((q_dot).^2); %%aceleracao de Xm

Lym = -Upm*sin(q)-Vpm*cos(q);%%coeficiente de aceleracao de Ym
Ym_dotdot = Kym.*(q_dotdot)+Lym.*((q_dot).^2); %%aceleracao de Ym

%%Biela
Lxb = -R*cos(q) + (-Upb*sin(A)+Vpb*cos(A)).*La + (Ka.^2).*(-Upb*cos(A) - Vpb*sin(A));%%coeficiente de aceleracao de Xb
Xb_dotdot = Kxb.*(q_dotdot)+Lxb.*((q_dot).^2); %%aceleracao de Xb

Lyb = -R*sin(q) + (-Upb*cos(A)-Vpb*sin(A)).*La + (Ka.^2).*(Upb*sin(A)-Vpb*cos(A));%%coeficiente de aceleracao de Xb
Yb_dotdot = Kyb.*(q_dotdot)+Lyb.*((q_dot).^2); %%aceleracao de Yb

%%Pistao
Yp_dotdot = 0; %%aceleracao de Yb
Xp_dotdot = X_dotdot; %%aceleracao de Xp


%%Animacao

%%Vetores posicao

Ypistao = (t.*Yp)./t; %%transformando em vetor Yp

P1 = [0;0]; %%ponto na junta fixa
P2 = [Xm;Ym]; %%ponto na junta de revolucao de ligacao biela-manivela
P3 = [Xp;Ypistao]; %%ponto na junta de revolucao de ligacao biela-pistao

y_pontos_pistao = [(11/117)*L (11/117)*L -(11/117)*L -(11/117)*L];

radii = 0.15;
% y = figure;
for i=1:length(t)
    
    ani = subplot(2,2,1);
    P1_ponto = viscircles([P1(1) P1(2)],radii);
    P2_ponto = viscircles([P2(1,i) P2(2,i)],radii);
    P3_ponto = viscircles([P3(1,i) P3(2,i)],radii);
    
    manivela_barra = line([P1(1) P2(1,i)],[P1(2) P2(2,i)],'LineWidth',2.5);
    biela_barra = line([P2(1,i) P3(1,i)],[P2(2,i) P3(2,i)],'LineWidth',2.5);
    
    x_pontos_pistao = Xp(i)+[(-28/117)*L (28/117)*L (28/117)*L (-28/117)*L];

    axis(ani,'equal');
    set(gca,'XLim',[-(120/117)*L (240/117)*L],'YLim',[-(120/117)*L (120/117)*L]);
    str = ['Time:' num2str(t(i)) 's'];
    tempo = text((220/117)*L,-(100/117)*L,str);
   
    hold on
    
    pistao = fill(x_pontos_pistao,y_pontos_pistao,'r');
    
    hold off
    
    view([90 -90])
    ylabel('Axis Y');
    xlabel('Axis X')
    
    pause(0.005);
    
    if i<length(t)
%         
        pos = subplot(2,2,2)
        yyaxis left
        plot(pos,t(1:i),Xp(1:i)) %%posicao X da junta biela-cursor
        xlabel('Time (s)')
        ylabel('Position Xp')
        hold on
        yyaxis right
        plot(pos,t(1:i),Yp_plot(1:i),'o-','MarkerSize',0.1)
        ylabel('Position Yp')
        hold off
        grid on
        
        pos1 = subplot(2,2,3)
        yyaxis left
        plot(pos1,t(1:i),Xm(1:i))
        xlabel('Time (s)')
        ylabel('Position Xm')
        hold on
        yyaxis right
        plot(pos1,t(1:i),Ym(1:i),'o-','MarkerSize',0.1)
        ylabel('Position Ym')
        hold off
        grid on
    
        pos2 = subplot(2,2,4)
        yyaxis left
        plot(pos2,t(1:i),Xb(1:i))
%         ylim([0 120])
        xlabel('Time (s)')
        ylabel('Position Xb')
        hold on
        yyaxis right
        plot(pos2,t(1:i),Yb(1:i),'o-','MarkerSize',0.1)
        ylabel('Position Yb')
        hold off
        grid on
        
        
        delete(P1_ponto)
        delete(P2_ponto)
        delete(P3_ponto)
        delete(manivela_barra)
        delete(biela_barra)
        delete(pistao)
        delete(tempo)
%         
    end
    
%      saveas(y,['Figura' num2str(i)],'png');
          
    
end