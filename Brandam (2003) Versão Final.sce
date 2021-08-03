// modelo de hidrolise enzimática de amido Brandam et al (2003)

clear
clc
tic()
dispvector=[1:1:100]
format(4)

// valores dos parâmetros
// gelatinização, T < Tg
k_g1 = 5.7*10^31         // s^-1
E_g1 = 220.6             // kJ/mol

// gelatinização, T > Tg
k_g2 = 3.1*10^14         // s^-1
E_g2 = 108.3             // kJ/mol
K = 273.15               // converter para°C
T_g = 60 + K             // ºC

// constante dos açúcares
k_gl = 0.023             // kg/U s
k_mlt = 0.117            // kg/U s
k_dex = 0.317            // kg/U s
k_alfa_mal = 0.389       // kg/U s
k_beta_mal = 0.137       // kg/U s
k_gl_ = 2.9*10^-8        // kg/U s
k_mlt_ = 1.5*10^-8       // kg/U s
k_alfa_mal_ = 1.2*10^-7  // kg/U s
k_beta_mal_ = 8.4*10^-8  // kg/U s

// constante dos Gases
R = 8.314/1000           // J/mol K

// massas molares
MM_gli = 180.156         // g/mol
MM_mal = 342.3           // g/mol
MM_mlt = 504.437         // g/mol
MM_dex = 504.4           // g/mol

// condições Iniciais
gasto = 0                // gasto inicial = 0
x1_0 = 98.51             // concentração amido sólido inicial em g/L
x2_0 = 0.0               // concentração amido gelatinizado inicial em g/L
x3_0 = 0.0               // concentração dextrinas inicial em g/L
x4_0 = 0.0               // concentração glicose inicial em g/L
x5_0 = 0.0               // concentração maltose inicial em g/L
x6_0 = 0.0               // concentração maltotriose inicial em g/L
Ts = [37,1,1,1]

function Taxa_Gelatinizacao = Sist_1(x1)
    if (T < T_g) then
        r_g = x1*k_g1*exp(-(E_g1/(R*T)))
    else
        r_g = x1*k_g2*exp(-(E_g2/(R*T)))
    end
    Taxa_Gelatinizacao = r_g
endfunction

function Atividade = Sist_2_3(T)
    if T < (63+273.15) then
        a_a = -0.001270154057*T^3 + 1.235930164837*T^2 - 400.438346*T + 43203.7983384673
        a_b = 0.049*T -13.9
    end
    
    if T >= (63+273.15) then
        a_a = 0.005646191991*T^3 - 5.814417*T^2 + 1995.013302*T - 228070.0234
        a_b = -0.374*T + 128.3
    end
    
    if a_a < 0 then
        a_a = 0
    end
    if a_b < 0 then
        a_b = 0
    end
    if T <313 then
        a_a = 1
        a_b = 1
    end
    a_Tref = 220.28986
    b_Tref = 486.95652
    Atividade = [a_a*a_Tref,a_b*b_Tref]
endfunction

function Amido = Sist_4(x2,T)
    
    a_ = Sist_2_3(T)
    a_a = a_(1)
    a_b = a_(2)
    
    r_gl = k_gl*a_a*x2
    r_mal = k_alfa_mal*a_a*x2+k_beta_mal*a_b*x2
    r_mlt = k_mlt*a_a*x2
    r_dex = k_dex*a_a*x2
    
    Amido = [r_gl,r_mal,r_mlt,r_dex]
endfunction

function Dextrinas = Sist_5(x3,T)
    
    a_ = Sist_2_3(T)
    a_a = a_(1)
    a_b = a_(2)
    
    r_gl_ = k_gl_*a_a*x3
    r_mal_ = k_alfa_mal_*a_a*x3 + k_beta_mal_*a_b*x3
    r_mlt_ = k_mlt_*a_a*x3

    Dextrinas = [r_gl_,r_mal_,r_mlt_]
endfunction

function out = f(in)
    
    for k=1:length(in)
        if in(k)<0
            in(k)=0
        end
    end
    
    x1 = in(1)
    x2 = in(2)
    x3 = in(3)
    
    r_g = Sist_1(x1)
    r_amido = Sist_4(x2,T)
    r_dextrinas = Sist_5(x3,T)
    
    r_gl = r_amido(1)
    r_mal = r_amido(2)
    r_mlt = r_amido(3)
    r_dex = r_amido(4)
    
    r_gl_ = r_dextrinas(1)
    r_mal_ = r_dextrinas(2)
    r_mlt_ = r_dextrinas(3)
    
    dx1 = -r_g //dSs_dt
    dx2 = r_g - r_gl - r_mal - r_mlt - r_dex //dSg_dt
    dx3 = r_dex - r_gl_ - r_mal_ - r_mlt_ //dD_dt
    dx4 = r_gl + r_gl_ //dgl_dt
    dx5 = r_mal + r_mal_; //dmal_dt
    dx6 = r_mlt + r_mlt_; //dmlt_dt
    
    out=[dx1,dx2,dx3,dx4,dx5,dx6]
endfunction 

in=[x1_0,x2_0,x3_0,x4_0,x5_0,x6_0]

function Cp = calor_especifico(T)
    
    Cp = (-4*10^(-11)*(T-K)^5 + 1*10^(-8)*(T-K)^4 - 1*10^(-6)*(T-K)^3 + 1*10^(-4)*(T-K)^2 - 0.0033*(T-K) + 4.2198) //kJ/kg K
    
endfunction

function Perfil_Temp = temperatura(t, Condicao)
    
    t=t/60
    if Condicao ==6 then
        Ts = [37,50,65,76] // teste 1 malt s1
        // Ts = [37,50,70,76] //teste 2 malt S1
        // Ts = [50,63,63,76] //teste 3 malt S1
        boundaryTimes = [0,20,57.5,77.5,20000000] //teste 1 malt s1 [5,5,20]
        //boundaryTimes = [0,45,90,225,20000000] //teste 3 malt S1 [5,5,25]
        risingTimes = [5,5,20]
        deltas=4
    end

    for k=1:(deltas)
        if t>=boundaryTimes(k) && t <boundaryTimes(k+1)
            T = Ts(k)
            if k>1  && t<(boundaryTimes(k)+ risingTimes(k-1))
                T = Ts(k-1) + (t-boundaryTimes(k))/risingTimes(k-1)* (Ts(k)- Ts(k-1))
        end
    end
end
    
Perfil_Temp = T + K
endfunction

dt = 0.3 // passo de integração
tf = 110*60 // tempo final em s
limit_time = tf/60

for t=dt:dt:tf
    
    percent_exec=(100*(t)/tf) 
    
     if isempty(dispvector(dispvector==percent_exec)) then
    else
        clc
        disp("Simulação "+ string(percent_exec)+"% concluída ...")
end
    
    T = temperatura(t, 6)
    
    k1 = f(in)
    k2 = f(in + 0.5*dt*k1)
    k3 = f(in + 0.5*dt*k2)
    k4 = f(in + dt*k3)
    
    for i=1:length(in)
        in(i) = in(i) + (dt/6)*(k1(i) + 2*k2(i) + 2*k2(i) + 2*k3(i) + k4(i))
    end
    
    tt(t/dt) = t/60
    TT(t/dt) = T-K
    d_x1(t/dt) = in(1) // amido 
    d_x2(t/dt) = in(2) // amido gelatinizado
    d_x3(t/dt) = in(3) // dextrinas
    d_x4(t/dt) = in(4) // glicose
    d_x5(t/dt) = in(5) // maltose
    d_x6(t/dt) = in(6) // maltotriose
    cp_plot(t/dt) = calor_especifico(T)

gasto = (gasto+((T+K)*dt)*(calor_especifico(T)))
produzido_gli = (d_x4(t/dt)-x4_0)
produzido_mal = (d_x5(t/dt)-x5_0)
balanco = (d_x1(t/dt)+d_x2(t/dt)+d_x3(t/dt)+d_x4(t/dt)+d_x5(t/dt)+d_x6(t/dt))-(x1_0+x2_0+x3_0+x4_0+x5_0+x6_0)
tempo_zero = (4.2)*(Ts(1)-25)*30*60

//concentrações finais das espécies em mol/L
//desejados
mols_gli = d_x4(t/dt)/MM_gli
mols_mal = d_x5(t/dt)/MM_mal

//indesejados
mols_dex = d_x3(t/dt)/MM_dex
mols_mlt = d_x6(t/dt)/MM_mlt

//seletividade
selet = (mols_gli+mols_mal)/(mols_mlt+mols_dex)

end

scf()

plot(tt,d_x1,'g.')
plot(tt,d_x3,'c.')
plot(tt,d_x4,'r.')
plot(tt,d_x5,'b.')
plot(tt,d_x6,'y.')
xlabel("Tempo (min)")
ylabel("Concentração (g/L)")
g=gca()
g.data_bounds =  [0,0;limit_time,100]

xtitle('Concentração ao longo do tempo',['Tempo (min)'],['Concentração (g/L)'])

h=legend(["Amido Sólido", "Dextrinas", "Glicose", "Maltose", "Maltotriose"], pos=1)

g1 = newaxes()
set(g1, "filled", "off")
plot(tt,TT,'k.',"markerSize",2)
ylabel("Temperatura (ºC)")
g1=gca()
g1.data_bounds =  [0,30;limit_time,80]
g1.axes_visible(1) = "off"
g1.y_location = "right"

clc
format(10)

disp("Gastando = "+ string((tempo_zero+gasto)/((limit_time+30)*60))+" kJ/kg")
disp("Rendimento de glicose = "+ string(produzido_gli)+" g/L")
disp("Rendimento de maltose = "+ string(produzido_mal)+" g/L")
disp("Balanço = "+ string(balanco)+" g/L")
disp("Seletividade de "+ string(selet))
disp("Simulacão 100% concluída em " + string(int(toc())) + " s")

