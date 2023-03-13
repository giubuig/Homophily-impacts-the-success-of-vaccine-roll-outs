library(data.table)
library(dplyr)
library(readxl)


# contact matrix ####
getContactMatrix <- function(country){
  if (country == "Catalunya"){
    ageS <- c(4.5, 5.2, 5.6, 5.2, 5.1, 5.7, 6.0, 7.1, 8.5, 8.3, 7.4, 6.6, 5.8, 5.0, 4.6, 3.6)/100
    dataC.dt <- read_excel("data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx", sheet = "Spain", col_names = F)
    dataC <- as.matrix(dataC.dt)
    NTot <- 7727029
  }
  
  if (country == "Italy"){
    ageS <-  c(3.7, 4.3, 4.8, 4.8, 5.0, 5.2, 5.5, 5.9, 6.8, 7.9, 8.1, 7.8, 6.7, 5.9, 5.8, 4.3)/100
    dataC.dt <- read_excel("data/contact_matrices_152_countries/MUestimates_all_locations_1.xlsx", sheet = "Italy", col_names = T)
    dataC <- as.matrix(dataC.dt)
    NTot <- 59257566
  }
  
  if (country == "Switzerland"){
    dataPop.dt <- read_excel("data/vaccination/Svizzera/svizzera_ageStructure.xlsx", col_names = F)
    dataPop.dt <- dataPop.dt[-c(1,2,84:114),1:2]
    dataPop.dt$...2 <- as.numeric(dataPop.dt$...2)
    NTot <- 8212011
    ageS <- rep(0,16)
    for (i in 1:16){
      ageS[i] <- sum(as.numeric(dataPop.dt$...2[(5*(i-1)+2):(i*5+1)]))/NTot
    }
    dataC.dt <- read_excel("data/contact_matrices_152_countries/MUestimates_all_locations_2.xlsx", sheet = "Switzerland", col_names = F)
    dataC <- as.matrix(dataC.dt)
  }
  
  if (country == "France"){
    dataPop.dt <- read_excel("data/vaccination/France/pop-totale-france.xls", col_names = F)
    dataPop.dt <- dataPop.dt[-c(1:6),c(2,5)]
    NTot <- 67063703
    ageS <- rep(0,16)
    for (i in 1:16){
      ageS[i] <- sum(as.numeric(dataPop.dt$...5[(5*(i-1)+1):(i*5)]))/NTot
    }
    dataC.dt <- read_excel("data/contact_matrices_152_countries/MUestimates_all_locations_1.xlsx", sheet = "France", col_names = T)
    dataC <- as.matrix(dataC.dt)
  }
  
  M <- 8
  ageN <- vector(length = M)
  for(i in 1:M){
    ageN[i] = ageS[2*(i-1)+1] + ageS[2*(i-1)+2]
  }
  N <- ageN * NTot
  
  dataCNew <- matrix(0.0, nrow = M, ncol = M)
  for(i in 1:M){
    for(j in 1:M){
      dataCNew[i,j] = (ageS[2*(i-1)+1]*(dataC[2*(i-1)+1,2*(j-1)+1] + dataC[2*(i-1)+1,2*(j-1)+2]) + 
                         ageS[2*(i-1)+2]*(dataC[2*(i-1)+2,2*(j-1)+1] + dataC[2*(i-1)+2,2*(j-1)+2]))/ageN[i] 
    }
  }
  
  #symmetrization
  C <- matrix(0.0, nrow = M, ncol = M)
  for(i in 1:M){
    for(j in 1:M){
      C[i,j] = (dataCNew[i,j]*N[i] + dataCNew[j,i]*N[j])/(2*N[i])
    }
  }
  
  return(list(N = N, C = C, NTot = NTot))
}

# vaccination data ####
### Catalunya
casesCat.dt <- fread("data/vaccination/Catalunya/catalunya_diari_total_pob.csv",header=TRUE)
casesCat.dt <- casesCat.dt[casesCat.dt$SEXE != "Altres",]
casesCat.dt <- casesCat.dt[casesCat.dt$GRUP_EDAT != "80 o m\xe9s",]
vacCat.dt <- as.data.table(casesCat.dt  %>%
                             dplyr::group_by(GRUP_EDAT, DATA) %>%
                             dplyr::summarise(vac = sum(VACUNATS_DOSI_2)))

colnames(vacCat.dt)[1:2] <- c("group","date")

totvacCat.dt <- as.data.table(vacCat.dt  %>%
                             dplyr::group_by(group) %>%
                             dplyr::summarise(tot_vac = sum(vac)))

### Italy
vacIta.dt <- fread("data/vaccination/Italia/somministrazioni-vaccini-latest.csv",header=TRUE)
vacIta.dt <- vacIta.dt[fascia_anagrafica != "80-89" & fascia_anagrafica != "90+"]
vacIta.dt <- vacIta.dt  %>%
  dplyr::group_by(fascia_anagrafica, data_somministrazione) %>%
  dplyr::summarise(vac = sum(seconda_dose))

l <- length(unique(vacIta.dt$data_somministrazione))
vacIta.dt <- rbind(data.table(fascia_anagrafica = rep("0-9",l), data_somministrazione = unique(vacIta.dt$data_somministrazione), vac = rep(0,l)), vacIta.dt)

colnames(vacIta.dt)[1:2] <- c("group","date")

totvacIta.dt <- as.data.table(vacIta.dt  %>%
                                dplyr::group_by(group) %>%
                                dplyr::summarise(tot_vac = sum(vac)))

### Switzerland
vacSwi.dt <- fread("data/vaccination/Svizzera/COVID19VaccPersons_AKL10_w_v2.csv",header=TRUE)
vacSwi.dt <- vacSwi.dt[altersklasse_covid19 != "80+" & altersklasse_covid19 != "65+" &
                         altersklasse_covid19 != "12 - 15" & altersklasse_covid19 != "16 - 64"]
vacSwi.dt <- vacSwi.dt[geoRegion == "CH" & type == "COVID19FullyVaccPersons"]
vacSwi.dt <- as.data.table(vacSwi.dt  %>%
  dplyr::group_by(altersklasse_covid19, date) %>%
  dplyr::summarise(vac = sum(entries)))

vacSwi.dt$date <- rep(seq(as.Date("2020-12-17"), as.Date("2021-11-10"), "weeks"),8)

colnames(vacSwi.dt)[1:2] <- c("group","date")

totvacSwi.dt <- as.data.table(vacSwi.dt  %>%
                                dplyr::group_by(group) %>%
                                dplyr::summarise(tot_vac = sum(vac)))

### France
vacFra.dt <- fread("data/vaccination/France/france.csv",header=TRUE)
vacFra.dt <- vacFra.dt[clage_vacsi != 0 & clage_vacsi != 80]
vacFra.dt <- as.data.table(vacFra.dt  %>%
                             dplyr::group_by(clage_vacsi, jour) %>%
                             dplyr::summarise(vac = sum(n_complet)))

l <- length(unique(vacFra.dt$jour))
vacFra.dt_new <- data.table(group = c(rep(9,l),rep(19,l),rep(29,l),rep(39,l),rep(49,l),rep(59,l),rep(69,l),rep(79,l)),
                            date = rep(seq(as.Date("2020-12-27"), as.Date("2021-10-06"), "days"), 8),
                            vac = rep(0,8*l))

vacFra.dt_new[group == 9]$vac <- vacFra.dt[clage_vacsi == 4]$vac + vacFra.dt[clage_vacsi == 9]$vac
vacFra.dt_new[group == 19]$vac <- vacFra.dt[clage_vacsi == 11]$vac + vacFra.dt[clage_vacsi == 17]$vac + round((2/7)*vacFra.dt[clage_vacsi == 24]$vac)
vacFra.dt_new[group == 29]$vac <- round((5/7)*vacFra.dt[clage_vacsi == 24]$vac) + vacFra.dt[clage_vacsi == 29]$vac
vacFra.dt_new[group == 39]$vac <- vacFra.dt[clage_vacsi == 39]$vac
vacFra.dt_new[group == 49]$vac <- vacFra.dt[clage_vacsi == 49]$vac
vacFra.dt_new[group == 59]$vac <- vacFra.dt[clage_vacsi == 59]$vac
vacFra.dt_new[group == 69]$vac <- vacFra.dt[clage_vacsi == 64]$vac + vacFra.dt[clage_vacsi == 69]$vac
vacFra.dt_new[group == 79]$vac <- vacFra.dt[clage_vacsi == 74]$vac + vacFra.dt[clage_vacsi == 79]$vac

totvacFra.dt <- as.data.table(vacFra.dt_new  %>%
                                dplyr::group_by(group) %>%
                                dplyr::summarise(tot_vac = sum(vac)))



# get homophily ####
getHomophily <- function(vac.dt, d0, d, N, C){

  new_vac.dt <- vac.dt[date >= d0 & date <= d]
  totvac.dt <- as.data.table(new_vac.dt  %>%
                               dplyr::group_by(group) %>%
                               dplyr::summarise(tot_vac = sum(vac)))
  V <- sum(totvac.dt$tot_vac)/sum(N)
  V_i <- totvac.dt$tot_vac/N
  
  N_contacts <- 0
  V_to_V <- 0
  N_to_N <- 0
  M <- length(V_i)
  for (i in 1:M){
    N_contacts <- N_contacts + N[i]*C[i,i]/2
    V_to_V <- V_to_V + V_i[i]*N[i]*C[i,i]*(V_i[i]-1/N[i])/2
    N_to_N <- N_to_N + (1-V_i[i])*N[i]*C[i,i]*(1-V_i[i]-1/N[i])/2
    j <- i + 1
    while (j <= M){
      N_contacts <- N_contacts + N[i]*C[i,j]
      V_to_V <- V_to_V + V_i[i]*N[i]*C[i,j]*V_i[j]
      N_to_N <- N_to_N + (1-V_i[i])*N[i]*C[i,j]*(1-V_i[j])
      j <- j + 1
    }
  }
  
  h = (V_to_V + N_to_N)/N_contacts
  a = (1 - h)/(2*V*(1-V))     #h = 1 - 2*a*V*(1-V)
  
  return(c(a, V))
}

getEvoHomophily <- function(vac.dt, d0, dates, N, C){
  evoHom.dt <- data.table(last_date = dates, homophily = rep(0,length(dates)), V = rep(0,length(dates)))
  
  for (i in 1:length(dates)){
    d <- dates[i]
    x <- getHomophily(vac.dt, d0, d, N, C)
    evoHom.dt[i,2] <- x[1]
    evoHom.dt[i,3] <- x[2]
  }
  
  return(evoHom.dt)
}


# results ####
### Catalunya
N_ageGroups <- getContactMatrix("Catalunya")[[1]]
K <- getContactMatrix("Catalunya")[[2]]

d0 <- min(vacCat.dt[vac>0]$date)   #day of first second dosis overall
dates <- unique(vacCat.dt[date >= d0 & date <= "2021-10-02"]$date)
evoHomCat.dt <- getEvoHomophily(vacCat.dt, d0, dates, N_ageGroups, K)

evoHomCat.dt_7d <- data.table(date = evoHomCat.dt$last_date[7:length(evoHomCat.dt$last_date)], h = 0, V = 0)
for(i in 1:dim(evoHomCat.dt_7d)[1]){
  evoHomCat.dt_7d$h[i] <- mean(evoHomCat.dt$homophily[i:(i+6)])
  evoHomCat.dt_7d$V[i] <- mean(evoHomCat.dt$V[i:(i+6)])
}

### Italy
N_ageGroups <- getContactMatrix("Italy")[[1]]
K <- getContactMatrix("Italy")[[2]]

d0 <- min(vacIta.dt[vac>0]$date)   #day of first second dosis overall
dates <- unique(vacIta.dt[date >= d0]$date)
evoHomIta.dt <- getEvoHomophily(vacIta.dt, d0, dates, N_ageGroups, K)

evoHomIta.dt_7d <- data.table(date = evoHomIta.dt$last_date[7:length(evoHomIta.dt$last_date)], h = 0, V = 0)
for(i in 1:dim(evoHomIta.dt_7d)[1]){
  evoHomIta.dt_7d$h[i] <- mean(evoHomIta.dt$homophily[i:(i+6)])
  evoHomIta.dt_7d$V[i] <- mean(evoHomIta.dt$V[i:(i+6)])
}

### Switzerland
N_ageGroups <- getContactMatrix("Switzerland")[[1]]
K <- getContactMatrix("Switzerland")[[2]]

d0 <- min(vacSwi.dt[vac>0]$date)   #day of first second dosis overall
dates <- unique(vacSwi.dt[date >= d0 & date <= "2021-10-02"]$date)
evoHomSwi.dt <- getEvoHomophily(vacSwi.dt, d0, dates, N_ageGroups, K)

### France
N_ageGroups <- getContactMatrix("France")[[1]]
K <- getContactMatrix("France")[[2]]

d0 <- min(vacFra.dt_new[vac>0]$date)   #day of first second dosis overall
dates <- unique(vacFra.dt_new[date >= d0 & date <= "2021-10-02"]$date)
evoHomFra.dt <- getEvoHomophily(vacFra.dt_new, d0, dates, N_ageGroups, K)

evoHomFra.dt_7d <- data.table(date = evoHomFra.dt$last_date[7:length(evoHomFra.dt$last_date)], h = 0, V = 0)
for(i in 1:dim(evoHomFra.dt_7d)[1]){
  evoHomFra.dt_7d$h[i] <- mean(evoHomFra.dt$homophily[i:(i+6)])
  evoHomFra.dt_7d$V[i] <- mean(evoHomFra.dt$V[i:(i+6)])
}
