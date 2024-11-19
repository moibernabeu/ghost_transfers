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

library(tidyverse)
library(patchwork)
library(see)

theme_set(theme_bw())
theme_replace(plot.tag = element_text(face = 'bold'))

tree_cols <- okabeito_colors(1:2)
names(tree_cols) <- c('ATPase', 'LUCA')

load_data <- function(x) {
  table <- shiftl[[x]]$table
  table$sim <- x
  return(table)
}

source('02_simulate_transfers.R')

# Stats about the branch space difference between both trees ----

# Branches in FECA-LECA
ATP_brsp <- read.csv('../outputs_atpase/node_ages_phylum.tsv', sep = '\t')
feca <- ATP_brsp[1, 'birth']
leca <- ATP_brsp[1, 'death']
ATP_brsp <- ATP_brsp[c(1,  which(ATP_brsp$birth >= leca & ATP_brsp$death <= feca & ATP_brsp$clades != 'Eukaryota')), ]
ATP_brspp <- ptr(ATP_brsp, segment_colour = tree_cols['ATPase'])

LUCA_brsp <- read.csv('../outputs_LUCA/node_ages_phylum.tsv', sep = '\t')
feca <- LUCA_brsp[1, 'birth']
leca <- LUCA_brsp[1, 'death']
LUCA_brsp <- LUCA_brsp[c(1,  which(LUCA_brsp$birth >= leca & LUCA_brsp$death <= feca & LUCA_brsp$clades != 'Eukaryota')), ]
LUCA_brsp[, 'birth'] = LUCA_brsp[, 'birth'] * 1000
LUCA_brsp[, 'death'] = LUCA_brsp[, 'death'] * 1000
LUCA_brspp <- ptr(LUCA_brsp, segment_colour = tree_cols['LUCA'])

dim(LUCA_brsp)[1] / dim(ATP_brsp)[1] * 100

mean(ATP_brsp[, 'length']) / mean(LUCA_brsp[, 'length'] * 1000)

# Age of the sister groups
sist_births <- rbind(data.frame(LUCA_brsp[-1, ], tree = 'LUCA'),
                     data.frame(ATP_brsp[-1, ], tree = 'ATPase'))
ggplot(sist_births, aes(-birth, colour = tree)) +
  stat_density(geom = 'line', position = 'identity') +
  geom_vline(xintercept = c(-feca * 1000, -leca * 1000), lty = 4, colour = 'steelblue') +
  xlab('Clades birth (Mya)') +
  ylab('Density') +
  scale_colour_manual(values = tree_cols) +
  coord_cartesian(xlim = c(-4500, 0)) +
  theme(legend.position = 'inside', legend.position.inside = c(0.9, 0.85))

# Loading data
load("../outputs_LUCA/shift_list.RData")
data_LUCA <- lapply(1:length(shiftl), load_data) %>% bind_rows()

load("../outputs_atpase/shift_list.RData")
data_atpase <- lapply(1:length(shiftl), load_data) %>% bind_rows()

data_LUCA$gt_diff <- (data_LUCA$ghost_birth / data_LUCA$transfer)
data_LUCA$gs_diff <- (data_LUCA$sister_birth / data_LUCA$ghost_birth)
data_LUCA$tree <- "LUCA"

data_atpase$gt_diff <- (data_atpase$ghost_birth / data_atpase$transfer)
data_atpase$gs_diff <- (data_atpase$sister_birth / data_atpase$ghost_birth)
data_atpase$tree <- "ATPase"

birth_transfer_diff <- rbind(data_LUCA, data_atpase)

birth_transfer_diff <- birth_transfer_diff %>%
  mutate(ghost_birth_tr = ifelse(tree == 'LUCA', ghost_birth * 1000, ghost_birth))

# Calculating the densities for the ghost births
LUCA_hist <- density(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'LUCA')],
                  from = min(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'LUCA')]), to = 4000) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = ifelse(x <= LUCA_brsp[1, 'birth'], 'Inside FECA-LECA', 'Outside FECA-LECA'),
         tree = 'LUCA')

ATP_hist <- density(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'ATPase')],
                     from = min(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'ATPase')]), to = 4000) %$% 
  data.frame(x = x, y = y) %>% 
  mutate(area = ifelse(x <= ATP_brsp[1, 'birth'], 'Inside FECA-LECA', 'Outside FECA-LECA'),
         tree = 'ATPase')

# Calculating the probability of a birth happening before FECA in both datasets
LUCA_ecdf <- ecdf(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'LUCA')])
ATP_ecdf <- ecdf(birth_transfer_diff$ghost_birth_tr[which(birth_transfer_diff$tree == 'ATPase')])

1 - LUCA_ecdf(LUCA_brsp[1, 'birth'])
1 - ATP_ecdf(ATP_brsp[1, 'birth'])

# Plotting the densities
c <- ggplot() +
  geom_ribbon(aes(-x, ymax = y, ymin = 0, fill = area, colour = tree, alpha = area), data = LUCA_hist) +
  geom_ribbon(aes(-x, ymax = y, ymin = 0, fill = area, colour = tree, alpha = area), data = ATP_hist) +
  geom_vline(data = data.frame(x = c(LUCA_brsp[1, 'birth'], LUCA_brsp[1, 'death'],
                                     ATP_brsp[1, 'birth'], ATP_brsp[1, 'death']),
                               tree = c('LUCA', 'LUCA', 'ATPase', 'ATPase')),
             aes(xintercept = -x, colour = tree), lty = 4) +
  scale_colour_manual(values = c(tree_cols)) +
  scale_alpha_manual(values = c('Inside FECA-LECA' = 0, 'Outside FECA-LECA' = 0.5)) +
  scale_fill_manual(values = c('Inside FECA-LECA' = NULL, 'Outside FECA-LECA' = 'grey80')) +
  xlab('Ghost birth') +
  ylab('Density') +
  coord_cartesian(xlim = c(-4000, -1500)) +
  theme(legend.position = 'none')

# Plotting the densities of the of the difference between the birth of the
# ghost and the transfer
d <- ggplot(birth_transfer_diff, aes(x = gt_diff, colour = tree)) +
  stat_density(position = 'identity', geom = 'line') +
  scale_colour_manual(values = tree_cols) +
  xlab('Ghost birth / transfer') +
  ylab('Density') +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.75, 0.7),
        legend.title = element_blank())

# Plotting the ghost birth vs detected donor ratio
ggplot(birth_transfer_diff, aes(gs_diff, colour = tree)) +
  stat_density(position = 'identity', geom = 'line') +
  scale_colour_manual(values = tree_cols) +
  xlab('Detected donor birth / ghost birth') +
  ylab('Density') +
  theme(legend.position = 'inside', legend.position.inside = c(0.75, 0.7))

# Generating the SF3 plot
pdf('../paper/SF_compbrl.pdf', width = 7, height = 7.75)
ATP_brspp + theme(axis.text.y = element_blank()) +
  LUCA_brspp +
  theme(axis.text.y = element_blank()) +
  c + theme(legend.position = 'none') +
  d +
  plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')') +
  plot_layout(heights = c(2.65, 1)) & theme(plot.tag = element_text(face = 'bold'))
dev.off()

