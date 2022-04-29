%% Topicos de Dinamica das Maquinas
%% Atividade

clc
close all
clear all

%%Declaracao dos parametros importantes
L = 117.2; %%mm
R = 40.1; %%mm
C = 0; %%mm

%%Massas e momentos de inercia
Mp = 0.2497; %%massa do pistao para aluminio (kg)
Mb = 0.12320; %%massa da biela (kg)
Mm = 1.725; %%massa da manivela (kg)
Ib = 526.079; %%momento de inercia da biela (kg.mm2)
Im = 7130.489; %%momento de inercia da manivela (kg.mm2) 

g = 9.81*(1000); %%mm/s2

% t = 0:0.0244:1; %%incremento do tempo (s)

t = 0:0.0244:(sqrt(0.4*pi)); %%incremente do tempo para MUV

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

q_dotdot = 10; %%aceleracao angular da generalizada (rad/s^2)
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

Upm_CM = 0; %%projecao do ponto de interesse em U (mm)
Vpm_CM = -1.016093; %%projecao do ponto de interesse em V (mm)

Xm = Upm*cos(q)-Vpm*sin(q); %%coordenada do ponto de interesse em X
Ym = Upm*sin(q)-Vpm*cos(q); %%coordenada do ponto de interesse em Y

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

Yp = C+Vp;
Xp = X+Up;

Vp_CM = 0.000332264; %%ponto de interesse em V (mm)
Up_CM = 13.1125; %%ponto de interesse em U (mm)

Yp_CM = C+Vp_CM;
Xp_CM = X+Up_CM;

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


%%Teorema dos eixos paralelos - Manivela

I_zero = Im + (Mm*(Upm^2 + Vpm^2)); %%momento de inercia da manivela

%%Equacao de Eksergian

Jq = I_zero + Mb*((Kxb.^2) + (Kyb.^2)) + Ib*(Ka.^2) + Mp*(Kx.^2); %%inercia generalizada 

Cq = (0.5)*( (Mb*(2*Kxb.*Lxb + 2*Kyb.*Lyb)) + 2*Ib*Ka.*La + 2*Mp*Kx.*Lx ); %%coeficiente centripeto

V = Mm*g*Xm + Mb*g*Xb + Mp*g*X; %%energia potencial

dVdq = Mm*g*Kxm + Mb*g*Kxb + Mp*g*Kx; %%derivada da energia potencial em relacao a q

%%Caso para massa do pistao diferente
Mp1 = 0.660; %%segundo caso para a massa do pistao (ferro fundido)

Cq1 = (0.5)*( (Mb*(2*Kxb.*Lxb + 2*Kyb.*Lyb)) + 2*Ib*Ka.*La + 2*Mp1*Kx.*Lx ); %%coeficiente centripeto
Jq1 = I_zero + Mb*((Kxb.^2) + (Kyb.^2)) + Ib*(Ka.^2) + Mp1*(Kx.^2); %%inercia generalizada 

V1 = Mm*g*Xm + Mb*g*Xb + Mp1*g*X; %%energia potencial

dVdq1 = Mm*g*Kxm + Mb*g*Kxb + Mp1*g*Kx; %%derivada da energia potencial em relacao a q


%%Energia Cinetica dos componentes

Tm = 0.5*(I_zero)*(q_dot.^2); %energia cinetica da manivela

Tm_plot = (t.*Tm)./t; %%vetorizacao de Tm

Tm_plot11 = 0.7*Tm_plot; %%vetorizacao de Tm

Tb = 0.5*(Mb*((Kxb.^2) + (Kyb.^2)) + Ib*(Ka.^2)).*(q_dot.^2); %energia cinetica da biela
Tp = 0.5*(Mp*(Kx.^2)).*(q_dot.^2); %energia cinetica do pistao

Tp1 = 0.5*(Mp1*(Kx.^2)).*(q_dot.^2); %energia cinetica do pistao (caso 2)

%%Energia potencial dos componentes

Vm = Mm*g*Xm; %%energia potencial da manivela
Vb = Mb*g*Xb; %%energia potencial da biela
Vp = Mp*g*X; %%energia potencial do pistao

Vp1 = Mp1*g*X; %%energia potencial do pistao (caso 2)



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

y_pontos_pistao = [11 11 -11 -11];

radii = 0.40;

y = figure;

for i=1:length(t)
    
    ani = subplot(2,2,1);
    P1_ponto = viscircles([P1(1) P1(2)],radii,'Color','k');
    P2_ponto = viscircles([P2(1,i) P2(2,i)],radii,'Color','r');
    P3_ponto = viscircles([P3(1,i) P3(2,i)],radii,'Color','k');
%     P4_ponto = viscircles([P4(1,i) P4(2,i)],radii,'Color','r');
    P5_ponto = viscircles([P5(1,i) P5(2,i)],radii,'Color','k');
    P6_ponto = viscircles([P6(1,i) P6(2,i)],radii);
    
    manivela_barra = line([P1(1) P3(1,i)],[P1(2) P3(2,i)],'LineWidth',2.5);
    biela_barra = line([P3(1,i) P5(1,i)],[P3(2,i) P5(2,i)],'LineWidth',2.5);
    
    P4_ponto = viscircles([P4(1,i) P4(2,i)],radii,'Color','r');
    
    x_pontos_pistao = Xp(i)+[(-28/117)*L (28/117)*L (28/117)*L (-28/117)*L];

    axis(ani,'equal');
    set(gca,'XLim',[-(120/117)*L (240/117)*L],'YLim',[-(120/117)*L (120/117)*L]);
    str = ['Lapsed time:' num2str(t(i)) 's'];
    tempo = text((220/117)*L,-(100/117)*L,str);
    
    
    hold on
    
    pistao = fill(x_pontos_pistao,y_pontos_pistao,[.7 .7 .7],'FaceAlpha',0.01);
    
    hold off
    
    view([90 -90])
    ylabel('Axis Y');
    xlabel('Axis X')
    
    pause(0.005);
    
    if i<length(t)
        
        pos = subplot(2,2,2)

        plot(q(1:i),Tm_plot11(1:i))
        hold on
        plot(q(1:i),Vm(1:i))
        hold off
        a = legend('Kinetic Energy - Crankshaft','Potential Energy - Crankshaft','Location','Best');
        a.FontSize = 7;
        ylabel('Energy (kg.mm^2/s^2)')
        xlabel('Angle q (rad)')
        ylim([-1*(10^8) 15*(10^8)])
        xlim([0 500])
%         ylim([-2*(10^5) 16*(10^5)])
%         xlim([0 2*pi])
        grid on

        
        pos1 = subplot(2,2,3)
        
        plot(q(1:i),Tb(1:i))
        hold on
        plot(q(1:i),Vb(1:i))
        hold off
        b = legend('Kinetic Energy - Rod','Potential Energy - Rod','Location','Best')
        b.FontSize = 7;
        ylabel('Energy (kg.mm^2/s^2)')
        xlabel('Angle q (rad)')
%         ylim([0 3*(10^7)])
%         xlim([0 500])
%         ylim([0 15*(10^4)])
%         xlim([0 2*pi])
        grid on
        
        
        pos2 = subplot(2,2,4)
        
        plot(q(1:i),Tp(1:i))
        hold on
%         plot(q(1:i),Tp1(1:i))
%         hold on
        plot(q(1:i),Vp(1:i))
        hold off
%         plot(q(1:i),Vp1(1:i))
%         hold off
        c = legend('Kinetic energy (Piston 1)','Potential energy (Piston 1)','Location','Best')
        c.FontSize = 7;
        ylabel('Energy (kg.mm^2/s^2)')
        xlabel('Angle q (rad)')
%         ylim([0 2.2*(10^8)])
%         xlim([0 500])
%         ylim([0 14*(10^5)])
%         xlim([0 2*pi])
        grid on
        
        
%         saveas(y,['Figura' num2str(i)],'png');
        
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
    
    
    

    
end




