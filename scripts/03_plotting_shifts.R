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
library(ggpubr)
library(tidyr)

theme_set(theme_bw())

load('../outputs/shift_list.RData')
shift <- gather(data.frame(do.call(rbind, lapply(shiftl, function(x) {x$shif_summary})), check.names = FALSE))
shift_older <- gather(data.frame(do.call(rbind, lapply(shiftl, function(x) {x$shift_and_older})), check.names = FALSE))

shift_older$key_2 <- factor(shift_older$key, labels = c(FALSE, TRUE))

a <- ggplot(shift, aes(x = key, y = value, colour = key, fill = key)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Shift') +
  ylim(10, 90)

b <- ggplot(shift_older, aes(x = key_2, y = value, colour = key, fill = key)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.2, alpha = 0) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Shift older than FECA') +
  ylim(10, 90)

pdf('../outputs/shift_older_proportions.pdf', width = 5, height = 4)
ggarrange(a, b, align = 'hv', legend = FALSE)
dev.off()

write.table(do.call(rbind, by(shift$value, shift$key, summary)),
            file = '../outputs/shift_summary.csv', sep = '\t')

load('../outputs/shift_list_ghosts_prop.RData')
shift <- data.frame(do.call(rbind, lapply(shiftl, function(x) {c(x$shift_summary, ghost_prop = x$ghosts_prop)})), check.names = FALSE) %>%
  pivot_longer(cols = c('No shift', 'Shift'))

ghosts_props_plot <- ggplot(shift, aes(x = as.factor(ghost_prop), y = value, colour = name, fill = name)) +
  geom_violin(alpha = 0.6, position = position_dodge(width = 1), draw_quantiles = c(0.5)) +
  # geom_boxplot(width = 0.2, alpha = 0, position = position_dodge(width = 1)) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab('Proportion') +
  xlab('Proportion of ghosts') +
  ylim(0, 100)

pdf('../outputs/shifts_proportion_ghosts.pdf', width = 6.75, height = 3)
ghosts_props_plot
dev.off()
