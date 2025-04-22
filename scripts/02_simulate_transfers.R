#!/usr/bin/env Rscript
#  Copyright Moisès Bernabeu, Saioa Manzano-Morales & Toni Gabaldón <moises.bernabeu.sci@gmail.com>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(ggplot2)
library(tidyr)
library(ggpubr)
library(stringr)
library(foreach)
library(doParallel)

theme_set(theme_bw())

ptr <- function(branch_space, sim_res = NULL, sim_res_2 = NULL, iter = NULL, segment_colour = 'black') {
  if (!is.null(sim_res_2)) {
    if (!is.null(iter)) {
      sim_res <- rbind(sim_res[iter, ], sim_res_2[iter, ])
    } else {
      sim_res <- rbind(sim_res, sim_res_2)
    }
  }
  
  feca <- branch_space[which(branch_space$node == 'FECA-LECA'), 'birth']
  leca <- branch_space[which(branch_space$node == 'FECA-LECA'), 'death']
  
  branch_space <- branch_space[-1, ]
  
  p <- ggplot(branch_space) +
    geom_segment(aes(x = -birth, xend = -death,
                     y = reorder(branch_space$node, birth),
                     yend = reorder(branch_space$node, birth)), colour = segment_colour) +
    geom_vline(xintercept = c(-feca, -leca),
               lty = 4, colour = 'steelblue') +
    xlab('Time (Mya)') +
    ylab('Branch')
  
    if (!is.null(sim_res)) {
      p <- p + geom_point(data = sim_res, aes(y = as.character(branch), x = -ghost_death), col = 'red') +
        geom_point(data = sim_res, aes(y = as.character(branch), x = -ghost_birth), col = 'steelblue') +
        geom_point(data = sim_res, aes(y = as.character(branch), x = -transfer), col = 'darkolivegreen')
    }
  
  print(p)
}

ghost_transfer <- function(branch_space, feca, leca) {
  # Random selection of a branch from the branch space
  br <- sample(1:dim(branch_space)[1], size = 1)
  donor <- branch_space[br, 'phylums']

  # Getting the selected branch birth and death
  birth <- branch_space[br, 'birth']
  death <- branch_space[br, 'death']

  # Defining the ghosts birth, transfer and death events
  if (birth <= feca & death >= leca) {
    # If the branch is inside FECA-LECA
    ghost_birth <- runif(1, death, birth)
    transfer <- runif(1, leca, ghost_birth)
    ghost_death <- runif(1, 0, transfer)
  } else if (birth >= feca & death >= leca) {
    # If the branch was born before FECA and died before LECA
    ghost_birth <- runif(1, death, birth)
    if (ghost_birth > feca) {
      transfer <- runif(1, leca, feca)
    } else {
      transfer <- runif(1, leca, ghost_birth)
    }
    ghost_death <- runif(1, 0, transfer)
  } else if (birth <= feca & death <= leca) {
    # If the branch was born after FECA and died after LECA
    ghost_birth <- runif(1, leca, birth)
    transfer <- runif(1, leca, ghost_birth)
    ghost_death <- runif(1, 0, transfer)
  } else if (birth >= feca & death <= leca) {
    # If the branch was born before FECA and died after LECA
    ghost_birth <- runif(1, leca, birth)
    if (ghost_birth > feca) {
      transfer <- runif(1, leca, feca)
    } else {
      transfer <- runif(1, leca, ghost_birth)
    }
    ghost_death <- runif(1, 0, transfer)
  }

  # Returning the information
  ghost <- data.frame('branch' = br,
                      'sister_birth' = birth, 'sister_death' = death,
                      'ghost_birth' = ghost_birth, 'ghost_death' = ghost_death,
                      'transfer' = transfer, 'donor' = donor)
  
  return(ghost)
}

ghost_simulator <- function(branch_space, feca, leca, N) {
  odf <- c()
  for (i in 1:N) {
    odf <- rbind(odf, ghost_transfer(branch_space, feca, leca))
  }
  return(odf)
}

get_shifts <- function(sim_1, sim_2, feca) {
  if (dim(sim_1)[1] != dim(sim_2)[1]) {
    stop('Both simulations have not the same number of iteration, they are not comparable.')
  }
  
  out <- c()
  for (iter in 1:1000) {
    pair <- rbind(sim_1[iter, ], sim_2[iter, ])
    mintrans <- which.min(pair$transfer)
    minbirth <- which.min(pair$ghost_birth)
    pair[mintrans, 'transfer_comparison'] <- 'early'
    pair[-mintrans, 'transfer_comparison'] <- 'late'
    pair[minbirth, 'conclusion'] <- 'early'
    pair[minbirth, 'shift'] <- pair[minbirth, 'transfer_comparison'] != pair[minbirth, 'conclusion']
    pair[-minbirth, 'conclusion'] <- 'late'
    pair[-minbirth, 'shift'] <- pair[-minbirth, 'transfer_comparison'] != pair[-minbirth, 'conclusion']
    pair[minbirth, 'older_th_feca'] <- pair[minbirth, 'ghost_birth'] > feca
    pair[-minbirth, 'older_th_feca'] <- pair[-minbirth, 'ghost_birth'] > feca
    pair[, 'iter'] <- iter

    row.names(pair) <- NULL
    out <- rbind(out, pair)
    # print(iter)
  }
  
  summ <- table(out$shift) / dim(out)[1] * 100
  names(summ) <- c('No shift', 'Shift')
  
  feca_summ <- table(out$older_th_feca) / dim(out)[1] * 100
  names(feca_summ) <- c('More recent than FECA', 'Older than FECA')

  shifts_id <- which(out$shift == TRUE)
  shift_and_older <- table(out[shifts_id, 'older_th_feca']) / length(shifts_id) * 100
  names(shift_and_older) <- c('Shift more recent than FECA', 'Shift older than FECA')

  return(list('table' = out, 'shif_summary' = summ, 'out_of_feca_summary' = feca_summ,
              'shift_and_older' = shift_and_older))
}

get_ghostprop_shifts <- function(sim_1, sim_2, feca, ghosts_prop = 1) {
  if (dim(sim_1)[1] != dim(sim_2)[1]) {
    stop('Both simulations have not the same number of iteration, they are not comparable.')
  }
  
  out <- c()
  for (iter in 1:dim(sim_1)[1]) {
    pair <- rbind(sim_1[iter, ], sim_2[iter, ])
    row.names(pair) <- NULL
    
    # Binomial to determine wether the detected ages are from the transfer (1)
    # or from the ghosth birth (0, meaning the absence of the closest relative)
    ghost <- rbinom(2, 1, 1 - ghosts_prop)
    
    ghost_birth_1 = pair[1, 'ghost_birth']
    transfer_1 = pair[1, 'transfer']
    donor_1 = pair[1, 'donor']
    branch_1 = pair[1, 'branch']
    if (ghost[1] == 1) {
      detected_1 = pair[1, 'transfer']
      older_1 = NA
    } else {
      detected_1 = pair[1, 'ghost_birth']
      older_1 = pair[1, 'ghost_birth'] > feca
    }
    
    ghost_birth_2 = pair[2, 'ghost_birth']
    transfer_2 = pair[2, 'transfer']
    donor_2 = pair[2, 'donor']
    branch_2 = pair[2, 'branch']
    if (ghost[2] == 1) {
      detected_2 = pair[2, 'transfer']
      older_2 = NA
    } else {
      detected_2 = pair[2, 'ghost_birth']
      older_2 = pair[2, 'ghost_birth'] > feca
    }
    
    shift <- (detected_1 >= detected_2) != (transfer_1 >= transfer_2)
    
    out <- rbind(out, data.frame(iter, branch_1, donor_1, transfer_1,
                                 ghost_birth_1, detected_1, older_1,
                                 branch_2, donor_2, transfer_2,
                                 ghost_birth_2, detected_2, older_2,
                                 shift))
  }
  
  summ <- table(out$shift) / dim(out)[1] * 100
  n <- c('FALSE' = 'No shift', 'TRUE' = 'Shift')
  names(summ) <- n[names(summ)]
  
  # Proportion of ghosts older than FECA
  olders <- na.omit(c(out$older_1, out$older_2))
  feca_summ <- table(olders) / length(olders) * 100
  names(feca_summ) <- c('More recent than FECA', 'Older than FECA')
  
  # Proportion of shifted conclusions with older than FECA ghost branch lengths
  olders <- na.omit(c(out[out$shift, 'older_1'], out[out$shift, 'older_2']))
  shift_and_older <- table(olders) / length(olders) * 100
  n <- c('FALSE' = 'Shift more recent than FECA', 'TRUE' = 'Shift older than FECA')
  names(shift_and_older) <- n[as.character(names(shift_and_older))]
  
  return(list('table' = out,
              'shift_summary' = summ,
              'out_of_feca_summary' = feca_summ,
              'shift_and_older' = shift_and_older,
              'ghosts_prop' = ghosts_prop))
}

# Loading data ----
brs <- read.csv('../outputs/node_ages_phylum.tsv', sep = '\t')

# Setting boundaries
feca <- brs[which(brs$node == 'FECA-LECA'), 'birth']
leca <- brs[which(brs$node == 'FECA-LECA'), 'death']

# Filtering donor branches
brs <- brs[which(brs$birth >= leca & brs$death <= feca & !str_detect(brs$clades, 'Eukaryota')), ]

registerDoParallel(112)
tm0 <- Sys.time()
shiftl <- foreach (i = 1:10000) %dopar% {
  get_shifts(ghost_simulator(brs, feca, leca, 1000), ghost_simulator(brs, feca, leca, 1000), feca)
}
tm1 <- Sys.time()
tm1 - tm0

save(shiftl, file = '../outputs/shift_list.RData')


registerDoParallel(4)
tm0 <- Sys.time()
shiftl_props <- c()
for (j in c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 1)) {
  propl <- foreach (i = 1:100) %dopar% {
    get_ghostprop_shifts(ghost_simulator(brs, feca, leca, 1000), ghost_simulator(brs, feca, leca, 1000), feca, j)
  }
  shiftl_props <- append(shiftl_props, propl)
}
tm1 <- Sys.time()
tm1 - tm0

save(shiftl_props, file = '../outputs/shift_list_ghosts_prop.RData')