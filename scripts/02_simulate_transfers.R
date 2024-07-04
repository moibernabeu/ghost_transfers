#!/usr/bin/env Rscript
#  Copyright Moisès Bernabeu, Saioa Manzano-Morales & Toni Gabaldón <saioa.manzano@bsc.es>
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

theme_set(theme_bw())

ptr <- function(branch_space, sim_res = NULL, sim_res_2 = NULL, iter = NULL) {
  if (!is.null(sim_res_2)) {
    if (!is.null(iter)) {
      sim_res <- rbind(sim_res[iter, ], sim_res_2[iter, ])
    } else {
      sim_res <- rbind(sim_res, sim_res_2)
    }
  }
  
  p <- ggplot(branch_space) +
    geom_segment(aes(x = -birth, xend = -death,
                     y = reorder(branch_space$node, birth),
                     yend = reorder(branch_space$node, birth))) +
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
    if (ghost_birth >= feca) {
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
  } else {
    # If the branch was born after FECA and died after LECA
    ghost_birth <- runif(1, leca, birth)
    transfer <- runif(1, leca, feca)
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

# Loading data ----
brs <- read.csv('../outputs/node_ages_phylum.tsv', sep = '\t')

# Setting boundaries
feca <- brs[which(brs$node == 'FECA-LECA'), 'birth']
leca <- brs[which(brs$node == 'FECA-LECA'), 'death']

# Filtering donor branches
brs <- brs[which(brs$birth >= leca & brs$death <= feca & brs$clades == 'Bacteria'), ]

# ptr(brs)

shiftl <- c()
tm0 <- Sys.time()
for (i in 1:10000) {
  shiftl[[i]] <- get_shifts(ghost_simulator(brs, feca, leca, 1000), ghost_simulator(brs, feca, leca, 1000), feca)
  print(i)
}
tm1 <- Sys.time()
tm1 - tm0

save(shiftl, file = '../outputs/shift_list.RData')
