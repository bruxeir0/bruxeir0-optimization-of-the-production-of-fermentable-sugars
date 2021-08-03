// código modelo de hidrolise enzimática de amido brandam et al (2003) anexado a otimização por enxame de partículas

clear
clc
tic()
format(10)

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
T_i = 62.562772  //PSO Foco em Gasto e Rendimento II
tempo_zero = (4.2)*(T_i-25)*30*60

// massas molares
MM_gli = 180.156         // g/mol
MM_mal = 342.3           // g/mol
MM_mlt = 504.437         // g/mol
MM_dex = 504.4           // g/mol

// condições iniciais
x1_0 = 98.51               //Concentração Amido Sólido Inicial em g/L
x2_0 = 0.0                 //Concentração Amido Gelatinizado Inicial em g/L
x3_0 = 0.0                 //Concentração Dextrinas Inicial em g/L
x4_0 = 0.0                 //Concentração Glicose Inicial em g/L
x5_0 = 0.0                 //Concentração Maltose Inicial em g/L
x6_0 = 0                 //Concentração Maltotriose Inicial em g/L

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

dt = 0.25 // passo de integração
tf = 60*60 // tempo final em s
limit_time = tf/60

function Cp = calor_especifico(T)
    
    Cp = (-4*10^(-11)*(T-K)^5 + 1*10^(-8)*(T-K)^4 - 1*10^(-6)*(T-K)^3 + 1*10^(-4)*(T-K)^2 - 0.0033*(T-K) + 4.2198) //kJ/kg K
    
endfunction

// PSO

function fitness = simulador(plottar,T0, c_ang,tempos)

    gasto = 0
    kant = 0
    b = T0
    T = T0 + 273.15
    t_min = 1e30
    gasto_min = 1e30
    tempos = [0,tempos,(tf-sum(tempos))]
    kk = length(tempos)
    tempos_original = tempos
    for k = 2:kk
        tempos(k)=sum(tempos_original(1:k))
    end

in=[x1_0,x2_0,x3_0,x4_0,x5_0,x6_0]

for t=dt:dt:tf-dt
    
    for k = 1:length(tempos)-1
    
        if t > tempos(k) && t <= (tempos(k+1))
            if k~=kant
                T0 = T - K
                b = T0 - c_ang(k)*t
            end
            T = c_ang(k)*t + b
            kant = k
        end
    end
    if T > T_max then
        T = T_max
    end
    if T < T_min then
        T = T_min
    end

    T=T+K
    
    k1 = f(in)
    k2 = f(in + 0.5*dt*k1)
    k3 = f(in + 0.5*dt*k2)
    k4 = f(in + dt*k3)
    
    for i=1:length(in)
        in(i) = in(i) + (dt/6)*(k1(i) + 2*k2(i) + 2*k2(i) + 2*k3(i) + k4(i))
    end
    
    tt(t/dt) = t/60
    TT(t/dt) = T-K
    d_x1(t/dt) = in(1) // Amido 
    d_x2(t/dt) = in(2) // Amido Gelatinizado
    d_x3(t/dt) = in(3) // Dextrinas
    d_x4(t/dt) = in(4) // Glicose
    d_x5(t/dt) = in(5) // Maltose
    d_x6(t/dt) = in(6) // Maltotriose
    
//gasto = gasto+(T*dt)
gasto = gasto+((T+K)*dt)*(calor_especifico(T))

evolucao_ami = (d_x1(t/dt))
produzido_gli = (d_x4(t/dt)-x4_0)
produzido_mal = (d_x5(t/dt)-x5_0)
produzido = (d_x4(t/dt)+d_x5(t/dt))-(x4_0+x5_0)

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



if plottar==1
    
disp("Gastando = "+ string((tempo_zero+gasto)/((limit_time+30)*60))+" kJ.s/kg")
disp("Gastando = "+ string(gasto)+" kJ.s/kg")
disp("Rendimento em glicose = "+ string(produzido_gli)+" g/L")
disp("Rendimento em maltose = "+ string(produzido_mal)+" g/L")
disp("Evolução do Amido Sólido = "+ string(evolucao_ami)+" g/L")
disp("Seletividade de "+ string(selet))
    
scf()
plot(tt,d_x1,'g.')
plot(tt,d_x3,'c.')
//plot(tt,d_x4,'r.')
//plot(tt,d_x5,'b.')
//plot(tt,d_x6,'y.',"markerSize",2)
xlabel("Time (min)")
ylabel("Concentration (g/L)")
g=gca()
g.data_bounds =  [0,0;limit_time,100]


xtitle('Concentração ao longo do tempo',['Time (min)'],['Concentration (g/L)'])

//h=legend(["Amido Sólido", "Dextrinas", "Glicose", "Maltose", "Maltotriose"], pos=-2)


g1 = newaxes()
set(g1, "filled", "off")
plot(tt,TT,'k.',"markerSize",2)
//xlabel("Tempo (min)")
ylabel("Temperature (ºC)")
g1=gca()
g1.data_bounds =  [0,25;limit_time,80]
g1.axes_visible(1) = "off"
g1.y_location = "right"
end

//if produzido < 78 then
//    produzido = 0
//end

    fitness = [produzido,gasto]
endfunction

tpop = 200
vezesmax = 60
c1 = 1
c2 = 2
T_min = 25
T_max = 350-K
T_gel = 63
T_a_max = 335-K
gasto_min =0
//gasto_max = (tempera_max*tf)

fim = 0

// partícula = [Tinicial. 5 c_angulares, 4 intervalos de tempo]

N_interv = 5
size_particula = 10
limiteinf = [T_min,-0.1,-0.1,-0.1,-0.1,-0.1,1,1,1,1]
limitesup = [T_max,0.1,0.1,0.1,0.1,0.1,tf,tf,tf,tf]
amp = limitesup - limiteinf
vezes = 1
inercia = 2
iner=ones(tpop,size_particula)*inercia

for i = 1:tpop
    for j = 1:size_particula
        velocidades(i,j) = rand()*amp(j)/10
        posicao(i,j)= rand()*amp(j) + limiteinf(j)
        if j > (N_interv + 2)
            tempousado = sum(posicao(i,7:(j-1)))
            tempolivre = tf - tempousado
            posicao(i,j) = rand()*tempolivre
        end
    end
    pbest(i,:) = [0,1e30]
    melhorposicao(i,:) = posicao(i,:)
end
gbest = [0,1e30]

posicao(1,:) = [64.668384  -0.0486105   0.0000459  -0.0011856   0.0385702   0.0190725   48.140101   3511.4108   8.5148174   0.2732467] //chute inicial

disp("entrou no loop principal")

while fim==0
    for i=1:tpop
        fitness(i,:) = simulador(0,posicao(i,1),posicao(i,2:6),posicao(i,7:10))
        
        if (fitness(i,1)> pbest(i,1))
            melhorposicao(i,:) = posicao(i,:)
            pbest(i,:) = fitness(i,:)
        end
        
        if (pbest(i,1) > gbest(1))
            gbest = pbest(i,:)
            melhorvez = vezes
            melhorpartic = i
            melhorposicaodetodos = posicao
        end
    printf('\n i=%d  vez=%d, particula= %.02f %.02f %.02f %.02f %.02f %.02f %.02f %.02f %.02f %.02f  \n sumtempo = %.01f, produzido = %.02f, gastando = %.02f \n Gbest = %.02f   %.02f  \n\n', i,vezes, posicao(i,1), posicao(i,2), posicao(i,3), posicao(i,4), posicao(i,5), posicao(i,6), posicao(i,7), posicao(i,8), posicao(i,9), posicao(i,10), sum(posicao(i,7:10)), fitness(i,1), fitness(i,2), gbest(1),gbest(2))
         
    end
    
    inercia=(inercia-inercia/2)
          if inercia<0.05 then
              inercia=2
          end
    
    for i=1:tpop
        for k=1:size_particula
            iner(i,k) = inercia
        velocidades(i,k)=inercia*velocidades(i,k)+ (c1* rand()*(melhorposicao(i,k)-posicao(i,k))) + (c2*rand()*melhorposicaodetodos(melhorpartic,k)-posicao(i,k))
        posicao(i,k)=posicao(i,k)+velocidades(i,k)
        if posicao(i,k)>limitesup(k) || posicao(i,k)<limiteinf(k)
            posicao(i,k)=rand()*amp(k)*0.3 + limiteinf(k)
        end
    end
        for k=1:size_particula

            if k>(N_interv+2) && sum(posicao(i,7:(k)))>=tf
            //disp(posicao(i,:))
            tempousado=sum(posicao(i,7:(k-1)))
            
            tempolivre=tf-tempousado
            if tempolivre<0
            disp(k)
            disp(tempousado)
            disp(tempolivre)
            disp(posicao(i,:))
            pause
            end    
            posicao(i,k)= rand()*tempolivre
            
            //disp(k)
            //disp(tempousado)
            //disp(tempolivre)
            //disp(posicao(i,:))
            //pause
        end
    end
    end
vetorgbest(vezes,:)=gbest
    vezes=vezes+1
    if vezes>=vezesmax then
        fim=1
    //end
end       

  if ((vezes/20)-int(vezes/20))<0.01 then
      
    for i = 1:tpop
        
    if i<>melhorpartic then        
        
    for j = 1:size_particula
        velocidades(i,j) = rand()*amp(j)/10
        posicao(i,j)= rand()*amp(j) + limiteinf(j)
        if j > (N_interv + 2)
            tempousado = sum(posicao(i,7:(j-1)))
            tempolivre = tf - tempousado
            posicao(i,j) = rand()*tempolivre
        end
    end
   end 
    //pbest(i,:) = [0,1e30]
    //melhorposicao(i,:) = posicao(i,:)
end
end
    
    
if (fim==1) then 
fitness(i,:) = simulador(1,melhorposicaodetodos(melhorpartic,1),melhorposicaodetodos(melhorpartic,2:6),melhorposicaodetodos(melhorpartic,7:10))
    
end 
end 

disp("Otimização concluída :)")
disp("")
disp("Tempo de execução: " + string(int(toc())) + "s")   
