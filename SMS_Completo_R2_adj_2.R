#C�digo SMS

setwd('C:/Users/fabri/Dropbox/UFF/2021/02 Segundo Semestre 2021/Disciplinas/Trabalho de Conclus�o de Curso II - Bacharelado/C�digos do R/SMS - Experimentos para o TCC/SVR_MSE_GA_R2_adj_2') #Define a pasta de arquivos onde ser�o salvos os arquivos de simula��o, de ordenamento pela RF, do corte pelo SVR e do refinamento pelo GA.

#install.packages('scrime') #instala o pacote SCRIME
library(scrime) #Carrega o pacote SCRIME
library(e1071)# Carrega o pacote e1071.
library(randomForest)
#library(doParallel)
#library(foreach)
library(GA) #Carrega a biblioteca do GA.
library(ggplot2) #Carrega a biblioteca ggplot2.

#Par�metros para simula��o

dados<-list();                   #Cria lista de dados.
genotipo<-list();                #Cria lista de gen�tipos.
num_individuos<-1000;            #N�mero de indiv�duos na amostra
num_snp<-100;                    #Total de marcadores simulados
list.snp<-list(c(1,2),c(3,4),c(5,6),c(7,8))  #Indica quais s�o os marcadores causais
list.ia<-list(c(1,2),c(3,-1),c(-2,-3),c(1,2))  #Contru��o das intera��es entre os 5 marcadores causais


#Executa a simula��o do gen�tipo e do fen�tipo simultaneamente
simulacao<-simulateSNPglm(n.obs=num_individuos,
                          n.snp=num_snp,
                          list.ia=list.ia,
                          list.snp=list.snp,
                          beta0=640,
                          beta=c(2,2,2,2),
                          maf=c(0.1,0.4),
                          err.fun=rnorm,
                          rand=123)

#Definindo o gen�tipo
genotipo[[1]]<-as.data.frame(simulacao$x)

#Definindo o fen�tipo
fenotipo<-as.data.frame(simulacao$y)
names(fenotipo)<-"fenotipo"
pdf(file="Histograma_fenotipo_2.pdf",height=5,width=9)
hist(fenotipo$fenotipo,xlab="Fen�tipo simulado",ylab="N�mero de touros",col="gray",main="")
dev.off()
pdf(file="Boxplot_fenotipo_2.pdf",height=5,width=9)
boxplot(fenotipo$fenotipo, ylab="Fen�tipo simulado")
dev.off()

#testando a normalidade do Fen�tipo
teste_normal<-shapiro.test(fenotipo$fenotipo) #N�o � normal

#definindo o dataframe genotipo com fenotipo
dados[[1]]<-as.data.frame(cbind(genotipo[[1]],fenotipo))
colnames(dados[[1]])[ncol(dados[[1]])]<-"fenotipo"

save(dados,file="Simulacao2.RData")

# Fun��o para realizar para o SVR com k-fold.
validacao_cruzada<-function(data,folds,gamma,cost,epsilon,kernel) 
{
  datalength = nrow(data)
  index <- 1:datalength
  size = trunc(datalength/folds)
  set.seed(123)
  geral = matrix(sample(index), ncol=folds, byrow=TRUE)
  
  mse = vector('double', folds)
  mape = vector('double', folds)
  corr = vector('double', folds)
  r2 = vector('double', folds)
  r2_adj = vector('double', folds)
  
  #  svm.pred<-foreach(i = 1:folds, 
  #                    .combine=rbind,
  #                    .inorder=TRUE,
  #                    .packages=c("e1071","rpart")) %dopar% { 
  
  for(i in 1:folds)
  {
    testv  = geral[,i];
    trainv = index[-testv]
    testset<-data.frame()
    trainset<-data.frame()
    
    testset  = as.data.frame(na.omit(data[testv,]))
    trainset = as.data.frame(na.omit(data[trainv,]))
    
    svmR.model <- svm(fenotipo~., data = trainset, kernel=kernel, gamma=gamma, cost=cost,epsilon=epsilon)
    #svmR.model <- svm(fenotipo~., data = trainset, kernel="linear", gamma=gamma, cost=cost,epsilon=epsilon)
    seltestset = as.data.frame(testset[,-ncol(testset)]);
    names(seltestset) = names(testset)[-ncol(testset)];
    svmR.pred <- predict(svmR.model, seltestset) 
    
    mse[i] = sum((svmR.pred-testset$fenotipo)^2)/dim(testset)[1]
    mape[i] = (sum( abs( (testset$fenotipo - svmR.pred) / testset$fenotipo) ) / length(testset$fenotipo)) * 100;
    corr[i]<-cor(svmR.pred,testset$fenotipo,method="pearson")
    r2[i] <- cor(svmR.pred,testset$fenotipo,method="pearson")^2
    r2_adj[i] <- 1-(1-r2[i])*((nrow(testset)-1)/(nrow(testset)-ncol(testset)-1))
  }
  
  
  
  resultado<-c(mean(mse),sd(mse),mean(mape),sd(mape),mean(corr),sd(corr),mean(r2),sd(r2),mean(r2_adj),sd(r2_adj))
  
  return(resultado)
}


########


#Fun��o para o c�lculo do valor p bruto e ajustado pela corre��o de Bonferroni
valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

################ In�cio do SMS##################

#Defini��o de listas para cada kernel do SVR
#svr_test<-list()
mean_svr_RF_list<-list() #Cria lista para a m�dia do SVR sobre o rank da RF
GA<-list()
minimo<-list()
corte<-list()
snps_selec_corte<-list()
snps_selec_ref<-list()
percentual_snps<-0.95;

i=1 #Contador do kernel do SVR utilizado

#Random Forest
ntree<-4000           #N�mero de �rvores dentro da floresta aleat�ria.

data_temp<-as.data.frame(dados[[1]]) #Transforma a base de dados em dataframe.

set.seed(1) #Semente aleat�ria para floresta aleat�ria.

RF<- randomForest(fenotipo~., data=data_temp,
                  ntree=ntree,
                  mtry=ncol(data_temp)-1,
                  importance=TRUE)    

rank_RF<-sort(importance(RF)[,1],decreasing=TRUE)
View(rank_RF)


#Erro MSE do SVR em rela��o ao rank da RF

#Inicializa��o de vari�vei do SVR
gamma = 0.01
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "linear"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#COnstr�i a sequ�ncia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequ�ncia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Gr�fico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_linear.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
rank_global_selec<-list()
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
View(rank_global)

#Selecionando somente o gen�tipo ap�s o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fen�tipo ap�s o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA sele��o dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[9];
  return(media);       
}
fitness<-function(x) {f(x)}


#Par�metros do GA
run = 30            #N�mero m�ximo de gera��es do GA sem melhoria.
maxiter = 10.000    #N�mero m�ximo de gera��es do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de muta��o.
elitism = 5         #N�mero de melhores indiv�duos do GA que ir�o para pr�xima gera��o sem altera��o alguma.
popSize=100         #Tamanho da popula��o em cada gera��o do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_linear_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_linear_GA_corr


#Relat�rio do AG
summary(GA[[i]])

#Gr�fico do AG
plot(GA[[i]])

#Salvando o gr�fico do AG com ggplot2
pdf(file="Grafico_GA_SVR_linear.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estat�sticas <- c(rep("Mediana", length(geracao)), rep("M�dia", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Gera��o = rep(geracao, 3),
  Aptid�o = c(mean_fitness, median_fitness, best_fitness),
  Estat�sticas = Estat�sticas
)

ggplot(data_grafico_1, aes(x = Gera��o, y = Aptid�o, group = Estat�sticas)) +
  geom_line(aes(colour = Estat�sticas, linetype = Estat�sticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Gera��o), max(data_grafico_1$Gera��o), by = 1))
dev.off()



############ SMS com SVR com kernel radial com gamma = 0.001 #############

i=2 #Contador do kernel "radial" com gamma=0.001 do SVR utilizado



#Erro MSE do SVR em rela��o ao rank da RF

#Inicializa��o de vari�vei do SVR
gamma = 0.001
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#COnstr�i a sequ�ncia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequ�ncia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Gr�fico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_radial_0001.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
View(rank_global)

#Selecionando somente o gen�tipo ap�s o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fen�tipo ap�s o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA sele��o dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[9];
  return(media);       
}
fitness<-function(x) {f(x)}


#Par�metros do GA
run = 30            #N�mero m�ximo de gera��es do GA sem melhoria.
maxiter = 10.000    #N�mero m�ximo de gera��es do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de muta��o.
elitism = 5         #N�mero de melhores indiv�duos do GA que ir�o para pr�xima gera��o sem altera��o alguma.
popSize=100         #Tamanho da popula��o em cada gera��o do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_radial_0001_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_radial_0001_GA_corr

#Relat�rio do AG
summary(GA[[i]])

#Gr�fico do AG
plot(GA[[i]])

#Salvando o gr�fico do AG com ggplot2
pdf(file="Grafico_GA_SVR_radial_0001.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estat�sticas <- c(rep("Mediana", length(geracao)), rep("M�dia", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Gera��o = rep(geracao, 3),
  Aptid�o = c(mean_fitness, median_fitness, best_fitness),
  Estat�sticas = Estat�sticas
)

ggplot(data_grafico_1, aes(x = Gera��o, y = Aptid�o, group = Estat�sticas)) +
  geom_line(aes(colour = Estat�sticas, linetype = Estat�sticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Gera��o), max(data_grafico_1$Gera��o), by = 1))
dev.off()





############ SMS com SVR com kernel radial com gamma = 0.01 #############

i=3 #Contador do kernel "radial" com gamma=0.01 do SVR utilizado

#Erro MSE do SVR em rela��o ao rank da RF

#Inicializa��o de vari�vei do SVR
gamma = 0.01
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#COnstr�i a sequ�ncia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequ�ncia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Gr�fico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_radial_001.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte,col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
View(rank_global)

#Selecionando somente o gen�tipo ap�s o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fen�tipo ap�s o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA sele��o dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[9];
  return(media);       
}
fitness<-function(x) {f(x)}


#Par�metros do GA
run = 30            #N�mero m�ximo de gera��es do GA sem melhoria.
maxiter = 10.000    #N�mero m�ximo de gera��es do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de muta��o.
elitism = 5         #N�mero de melhores indiv�duos do GA que ir�o para pr�xima gera��o sem altera��o alguma.
popSize=100         #Tamanho da popula��o em cada gera��o do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_radial_001_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_radial_001_GA_corr

#Relat�rio do AG
summary(GA[[i]])

#Gr�fico do AG
plot(GA[[i]])

#Salvando o gr�fico do AG com ggplot2
pdf(file="Grafico_GA_SVR_radial_001.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estat�sticas <- c(rep("Mediana", length(geracao)), rep("M�dia", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Gera��o = rep(geracao, 3),
  Aptid�o = c(mean_fitness, median_fitness, best_fitness),
  Estat�sticas = Estat�sticas
)

ggplot(data_grafico_1, aes(x = Gera��o, y = Aptid�o, group = Estat�sticas)) +
  geom_line(aes(colour = Estat�sticas, linetype = Estat�sticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Gera��o), max(data_grafico_1$Gera��o), by = 1))
dev.off()



############ SMS com SVR com kernel radial com gamma = 0.1 #############

i=4 #Contador do kernel "radial" com gamma=0.01 do SVR utilizado

#Erro MSE do SVR em rela��o ao rank da RF

#Inicializa��o de vari�vei do SVR
gamma = 0.1
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#COnstr�i a sequ�ncia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequ�ncia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Gr�fico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_radial_01.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte[[i]],col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
View(rank_global)

#Selecionando somente o gen�tipo ap�s o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fen�tipo ap�s o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA sele��o dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[9];
  return(media);       
}
fitness<-function(x) {f(x)}


#Par�metros do GA
run = 30            #N�mero m�ximo de gera��es do GA sem melhoria.
maxiter = 10.000    #N�mero m�ximo de gera��es do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de muta��o.
elitism = 5         #N�mero de melhores indiv�duos do GA que ir�o para pr�xima gera��o sem altera��o alguma.
popSize=100         #Tamanho da popula��o em cada gera��o do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_radial_01_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_radial_01_GA_corr

#Relat�rio do AG
summary(GA[[i]])

#Gr�fico do AG
plot(GA[[i]])

#Salvando o gr�fico do AG com ggplot2
pdf(file="Grafico_GA_SVR_radial_01.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estat�sticas <- c(rep("Mediana", length(geracao)), rep("M�dia", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Gera��o = rep(geracao, 3),
  Aptid�o = c(mean_fitness, median_fitness, best_fitness),
  Estat�sticas = Estat�sticas
)

ggplot(data_grafico_1, aes(x = Gera��o, y = Aptid�o, group = Estat�sticas)) +
  geom_line(aes(colour = Estat�sticas, linetype = Estat�sticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Gera��o), max(data_grafico_1$Gera��o), by = 1))
dev.off()


############ SMS com SVR com kernel radial com gamma = 1 #############

i=5 #Contador do kernel "radial" com gamma=0.01 do SVR utilizado

#Erro MSE do SVR em rela��o ao rank da RF

#Inicializa��o de vari�vei do SVR
gamma = 1
cost = 1.0
epsilon = 0.1
folds = 10
kernel = "radial"

#Etapa de Corte com o SVR sobre o Rank da RF
mean_svr_RF_list[[i]]<-vector()
passo<-10;

limite<-floor((length(names(dados[[1]][,-dim(dados[[1]])[2]]))/passo)*percentual_snps);


#COnstr�i a sequ�ncia crescente de SNPs mais importantes a partir do rank da RF
for (cont in 1:(limite+1))
{
  
  print(cont)
  j<-cont*passo
  #Ajuste quando j==1, pois o R demonstra um erro quando j==1.
  if (j==1) {
    
    var_sel<-names(rank_RF)[1]
    
  }else{
    
    var_sel<-names(rank_RF)[1:j]
    
  }
  
  #Avalia os 10-fold para o SVR na sequ�ncia crescente de SNPs a partir do rank da RF
  svr_test<- validacao_cruzada(data = dados[[1]][c(var_sel,"fenotipo")],
                               folds = folds, 
                               gamma = gamma,
                               cost = cost,
                               epsilon = epsilon, 
                               kernel = kernel)
  mean_svr_RF_list[[i]][cont]<-svr_test[1]
} 

#Constroi o Gr�fico do MSE do SVR sobre a RF
pdf(file="Grafico_MSE_SVR_radial_1.pdf",height=5,width=9)
plot(seq(passo,limite*(passo)+10,by=passo),
     mean_svr_RF_list[[i]],
     type="o",
     lwd=2,
     xlab="Grupo de Marcadores",
     ylab="MSE do SVR")

#ESCOLHE O PRIMEIRO PONTO DE CORTE
minimo[[i]]<-which.min(mean_svr_RF_list[[i]])
corte[[i]]<-(minimo[[i]]+1)*passo #Adotou-se o segundo menor MSE do SVR, pois o menor foi o referente ao primeiro marcador SNP2
abline(v=corte,col="black",lty=2)
dev.off()


snps_selec_corte[[i]]<-names(rank_RF[1:corte[[i]]])
snps_selec_corte[[i]]

#Construindo o dataframe com rank_RF juntamente com o p_valor dos SNPs selecionados na etapa de corte
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]][c(snps_selec_corte[[i]],'fenotipo')])
rank_global_selec[[i]]<-cbind(rank_RF[1:corte[[i]]],valor_p)
View(rank_global_selec[[i]])

#Construindo o dataframe com rank_RF juntamente com o p_valor de todos os SNPs da base de dados inicial.
valor_p<-data.frame();
valor_p<-valor.p(dados[[1]])
rank_global<-cbind(importance(RF)[,1],valor_p)
View(rank_global)

#Selecionando somente o gen�tipo ap�s o corte
genotipo[[2]]<-dados[[1]][,names(rank_RF[1:corte[[i]]])]

#Selecionando o genotipo e o fen�tipo ap�s o corte
dados[[2]]<-cbind(genotipo[[2]],dados[[1]]$fenotipo)
names(dados[[2]])[ncol(dados[[2]])]<-"fenotipo"


###AG para a SEGUNDA sele��o dos marcadores na base de dados inicial (REFINAMENTO)

#Declara o objeto GA como uma lista

f<-function(x)
{
  inc<-which(x==1)
  
  cat("inc = ",inc, "\n")
  
  dados_validacao<- dados[[2]][,c(inc,ncol(dados[[2]]))]
  model<-validacao_cruzada(dados_validacao,folds,gamma,cost,epsilon,kernel);
  media<-model[9];
  return(media);       
}
fitness<-function(x) {f(x)}


#Par�metros do GA
run = 30            #N�mero m�ximo de gera��es do GA sem melhoria.
maxiter = 10.000    #N�mero m�ximo de gera��es do GA.
pcross = 0.8        #Probabilidade de crossing over
pmut = 0.1          #Probabilidade de muta��o.
elitism = 5         #N�mero de melhores indiv�duos do GA que ir�o para pr�xima gera��o sem altera��o alguma.
popSize=100         #Tamanho da popula��o em cada gera��o do GA.

GA[[i]]<-ga(type="binary",
            fitness=fitness, 
            nBits=ncol(genotipo[[2]]),
            popSize=popSize,
            names=colnames(genotipo[[2]]),
            maxiter=maxiter,
            seed=i,
            parallel=TRUE,
            run=run, 
            suggestions=matrix(rep(1,ncol(genotipo[[2]])),
                               ncol=ncol(genotipo[[2]])))

snps_selec_ref[[i]]<-names(GA[[i]]@solution[,which(GA[[i]]@solution==1)])

snps_selec_SVR_radial_1_GA_corr<-snps_selec_ref[[i]]
snps_selec_SVR_radial_1_GA_corr

#Relat�rio do AG
summary(GA[[i]])

#Gr�fico do AG
plot(GA[[i]])

#Salvando o gr�fico do AG com ggplot2
pdf(file="Grafico_GA_SVR_radial_1.pdf",height=5,width=9)
geracao<-seq(1,GA[[1]]@iter,by=1)
mean_fitness   <-  GA[[i]]@summary[,2]
median_fitness <- GA[[i]]@summary[,4]
best_fitness   <- GA[[i]]@summary[,1]
Estat�sticas <- c(rep("Mediana", length(geracao)), rep("M�dia", length(geracao)),rep("Melhor", length(geracao)) )
data_grafico_1 <- data.frame(
  Gera��o = rep(geracao, 3),
  Aptid�o = c(mean_fitness, median_fitness, best_fitness),
  Estat�sticas = Estat�sticas
)

ggplot(data_grafico_1, aes(x = Gera��o, y = Aptid�o, group = Estat�sticas)) +
  geom_line(aes(colour = Estat�sticas, linetype = Estat�sticas),size=2) + 
  geom_point() +
  scale_x_continuous(breaks = seq(min(data_grafico_1$Gera��o), max(data_grafico_1$Gera��o), by = 1))
dev.off()

####Uni�o dos SNPs selecionados####
#O argumento lista_SNPs tem que est� no formato lista com pelo menos 2 elementos.
uniao_snps<-function(lista_SNPs){
  uniao <- list()
  uniao[[1]]<-union(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    uniao[[a+1]]<-union(uniao[[a]],lista_SNPs[[a+2]])
  }
  uniao_final <- uniao[[a+1]]
  return(uniao_final)
}

uniao_snps(snps_selec_ref)

uniao_final<-uniao_snps(snps_selec_ref)

####Interse��o dos SNPs selecionados####
#O argumento lista_SNPs tem que est� no formato lista com pelo menos 2 elementos.
intersecao_snps<-function(lista_SNPs){
  intersecao <- list()
  intersecao[[1]]<-intersect(lista_SNPs[[1]],lista_SNPs[[2]])
  for (a in 1:(length(lista_SNPs)-2)){
    intersecao[[a+1]]<-intersect(intersecao[[a]],lista_SNPs[[a+2]])
  }
  intersecao_final <- intersecao[[a+1]]
  if (length(intersecao_final) == 0L) {
    return("Conjunto vazio")
  } else
    return(intersecao_final)
}

intersecao_snps(snps_selec_ref)

intersecao_final<-intersecao_snps(snps_selec_ref)


#####Valor-p bruto e ajustado - Sele��o#####
#C�lculo do valor p bruto e ajustado pela corre��o de Bonferroni
valor.p<-function(genotipo_fenotipo)
{
  valor_p_bruto<-vector();
  valor_p_ajustado<-vector();
  y<-genotipo_fenotipo[,ncol(genotipo_fenotipo)];
  m<-ncol(genotipo_fenotipo)-1;
  model_regression<-list();
  saida<-data.frame();
  
  for (i in 1:(ncol(genotipo_fenotipo)-1))
  {
    x<-genotipo_fenotipo[,i];
    model_regression[[i]]<-lm(y~x,data=genotipo_fenotipo);
    valor_p_bruto[i]<-ifelse(all(x==3)|all(x==2)|all(x==1),1,summary(model_regression[[i]])[[4]][2,4]);
    
  }
  valor_p_ajustado<- m*valor_p_bruto;  
  valor_p<-data.frame();
  valor_p<-cbind(valor_p_bruto,valor_p_ajustado);
  colnames(valor_p)<-c("Valor p bruto","Valor p ajustado")
  rownames(valor_p)<-names(genotipo_fenotipo)[1:(ncol(genotipo_fenotipo)-1)];
  return(valor_p)
}

valor_p<-as.data.frame(valor.p(dados[[1]]))
valor_p_order<-valor_p[order(valor_p[,1], decreasing = FALSE),]
View(valor_p_order)
valor_p_bruto_order <- valor_p_order[valor_p_order$`Valor p bruto`<= 0.05,] # Filtrando somente os SNPs com valor-p bruto <= 0.05
valor_p_corrigido_order <- valor_p_order[valor_p_order$`Valor p ajustado`<= 0.05,] # Filtrando somente os SNPs com valor-p corrigido <= 0.05
valor_p_bruto_selecao<-rownames(valor_p_bruto_order) # SNPs selecionados pelo valor-p bruto <= 0.05
valor_p_corrigido_selecao<-rownames(valor_p_corrigido_order) # SNPs selecionados pelo valor-p bruto <= 0.05
valor_p_bruto_selecao # Exibe os SNPs selecionados pelo valor-p bruto <= 0.05
valor_p_corrigido_selecao # Exibe os SNPs selecionados pelo valor-p corrigido <= 0.05


#Imprimindo as informa��es finais para o texto

cat('SMS Linear =',snps_selec_ref[[1]])

cat('SMS Radial gamma 0,001 =',snps_selec_ref[[2]])

cat('SMS Radial gamma 0,01 =',snps_selec_ref[[3]])

cat('SMS Radial gamma 0,1 =',snps_selec_ref[[4]])

cat('SMS Radial gamma 1 =',snps_selec_ref[[5]])

cat('Uni�o =',uniao_final)

cat('Interse��o =',intersecao_final)

cat('Valor-p bruto =',valor_p_bruto_selecao)

cat('Valor-p corrigido =',valor_p_corrigido_selecao)

################# Medidas de Otimalidade dos SNPs causais###################

var_sel_causais<-c("SNP1","SNP2","SNP3","SNP4","SNP5","SNP6","SNP7","SNP8")

resultado_causais_linear<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                            folds = folds, 
                                            gamma = 0.001,
                                            cost = cost,
                                            epsilon = epsilon, 
                                            kernel = "linear")

resultado_causais_radial_0.001<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                  folds = folds, 
                                                  gamma = 0.001,
                                                  cost = cost,
                                                  epsilon = epsilon, 
                                                  kernel = "radial")

resultado_causais_radial_0.01<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                 folds = folds, 
                                                 gamma = 0.01,
                                                 cost = cost,
                                                 epsilon = epsilon, 
                                                 kernel = "radial")

resultado_causais_radial_0.1<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                                folds = folds, 
                                                gamma = 0.1,
                                                cost = cost,
                                                epsilon = epsilon, 
                                                kernel = "radial")

resultado_causais_radial_1<-validacao_cruzada(data = dados[[1]][c(var_sel_causais,"fenotipo")],
                                              folds = folds, 
                                              gamma = 1,
                                              cost = cost,
                                              epsilon = epsilon, 
                                              kernel = "radial")
###########################################################################
#Sa�da do SVR kernel linear com somente SNPs causais
cat('Correla��o do SVR Linear = ',resultado_causais_linear[5])
cat('R2 ajustado do SVR Linear = ',resultado_causais_linear[9])
cat('MSE do SVR Linear = ',resultado_causais_linear[1])

#Sa�da do SVR kernel radia gamma = 0.001 com somente SNPs causais
cat('Correla��o do SVR Radial 0,001 = ',resultado_causais_radial_0.001[5])
cat('R2 ajustado do SVR Radial 0,001 = ',resultado_causais_radial_0.001[9])
cat('MSE do SVR Radial 0,001 = ',resultado_causais_radial_0.001[1])

#Sa�da do SVR kernel radia gamma = 0.01 com somente SNPs causais
cat('Correla��o do SVR Radial 0,01 = ',resultado_causais_radial_0.01[5])
cat('R2 ajustado do SVR Radial 0,01 = ',resultado_causais_radial_0.01[9])
cat('MSE do SVR Radial 0,01 = ',resultado_causais_radial_0.01[1])

#Sa�da do SVR kernel radia gamma = 0.1 com somente SNPs causais
cat('Correla��o do SVR Radial 0,1 = ',resultado_causais_radial_0.1[5])
cat('R2 ajustado do SVR Radial 0,1 = ',resultado_causais_radial_0.1[9])
cat('MSE do SVR Radial 0,1 = ',resultado_causais_radial_0.1[1])

#Sa�da do SVR kernel radia gamma = 1 com somente SNPs causais
cat('Correla��o do SVR Radial 1 = ',resultado_causais_radial_1[5])
cat('R2 ajustado do SVR Radial 1 = ',resultado_causais_radial_1[9])
cat('MSE do SVR Radial 1 = ',resultado_causais_radial_1[1])

save.image(file="SMS_R2_adj_2.RData")
