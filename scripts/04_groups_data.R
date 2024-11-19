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

## Libraries import
library(tidyverse)
library(pheatmap)

## Data import
load('../outputs_LUCA/shift_list.RData')

f <- function(x, n, i) {
  print(i)
  return(data.frame(simulation = n[[i]], shiftl[[i]][['table']][, c('iter', 'donor', 'shift')]))
}

simulations <- do.call(rbind, lapply(seq_along(shiftl), f, x = shiftl, n = 1:length(shiftl)))

write.table(simulations, '../outputs_LUCA/all_simulations.tsv', sep = '\t', quote = FALSE)

donor_pairs <- simulations %>%
  group_by(simulation, iter) %>%
  mutate(donor_no = order(shift)) %>%
  pivot_wider(names_from = donor_no,
              values_from = donor,
              names_prefix = 'donor_')

real_donors <- unique(c(donor_pairs$donor_1, donor_pairs$donor_2))
if (sum(str_detect(real_donors, ';')) >= length(LETTERS)) {
  LETTERS <- c(paste0('A', LETTERS), paste0('B', LETTERS))
}
donors <- c()
anc = 1
for (i in 1:length(real_donors)) {
  if (str_detect(real_donors[i], ';')) {
    donors <- c(donors, paste0('Ancestor ', LETTERS[anc]))
    anc = anc + 1
  } else {
    donors <- c(donors, real_donors[i])
  }
}

names(donors) <- real_donors
donor_correspondences <- data.frame(donors, real_donors) %>%
  separate_longer_delim(real_donors, ';')
write.table(donor_correspondences, '../outputs_LUCA/ancestors_correspondence.tsv', sep = '\t', quote = FALSE)

total <- matrix(nrow = length(donors), ncol = length(donors), dimnames = list(donors, donors))
shift <- matrix(nrow = length(donors), ncol = length(donors), dimnames = list(donors, donors))
done <- c()
for (i in real_donors) {
  for (j in real_donors) {
    if (!paste0(j, i) %in% done) {
      done <- c(done, paste(i, j))
      if (i != j) {
        total_df <- donor_pairs[which((donor_pairs$donor_1 == i & donor_pairs$donor_2 == j) |
                                        (donor_pairs$donor_1 == j & donor_pairs$donor_2 == i)), ] %>%
          na.omit()
        shift[donors[i], donors[j]] <- sum(total_df[, 'shift'], na.rm = TRUE)
        shift[donors[j], donors[i]] <- sum(total_df[, 'shift'], na.rm = TRUE)
        total[donors[j], donors[i]] <- dim(total_df)[1]
        total[donors[j], donors[i]] <- dim(total_df)[1]
        
      } else {
        shift[donors[i], donors[j]] <- 0
        total[donors[i], donors[j]] <- 0
      }
    }
  }
}

write.table(shift / total, file = '../outputs_LUCA/shiftprop_clades_matrix.tsv',
            sep = '\t', row.names = TRUE, col.names = TRUE, quote = FALSE)

pheatmap(shift / total, na_col = 'white', cellheight = 12, cellwidth = 12,
         breaks = 0:100/100, cluster_cols = 2, cluster_rows = 2, cutree_rows = 4,
         cutree_cols = 4, filename = '../outputs_LUCA/heatmap.pdf')
