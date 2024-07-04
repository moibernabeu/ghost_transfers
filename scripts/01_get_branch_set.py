#!/usr/bin/env python3
'''
Copyright Moisès Bernabeu, Saioa Manzano-Morales & Toni Gabaldón <saioa.manzano@bsc.es>

This file is free software: you may copy, redistribute and/or modify it  
under the terms of the GNU General Public License as published by the  
Free Software Foundation, either version 2 of the License, or (at your  
option) any later version.  

This file is distributed in the hope that it will be useful, but  
WITHOUT ANY WARRANTY; without even the implied warranty of  
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
General Public License for more details.  

You should have received a copy of the GNU General Public License  
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

import ete3

treefile = '../data/pruned_tree.nwk'

t = ete3.PhyloTree(treefile)

euk_leaves = []
for leaf in t.iter_leaves():
    leaf.add_feature('clade', leaf.name.split('_', 1)[0].replace("'", ''))
    leaf.add_feature('phylum', leaf.name.split('_', 1)[1].split('-', 1)[0].split('_', 1)[0].replace("'", ''))
    if '-' in leaf.name:
        leaf.add_feature('classes', leaf.name.split('_', 1)[1].split('-')[1].split('_', 1)[0].replace("'", ''))
    else:
        leaf.add_feature('classes', leaf.name.split('_', 1)[1].split('_')[1].split('_', 1)[0].replace("'", ''))
    if leaf.clade == 'Eukaryota':
        euk_leaves.append(leaf.name)

leca = t.get_common_ancestor(euk_leaves)
leca_age = leca.get_distance(leca.get_leaf_names()[0])
laeca_age = leca_age + leca.dist

ofile = open('../outputs/node_ages_phylum.tsv', 'w')
ofile.write('\t'.join(['node', 'birth', 'death', 'length', 'clades', 'phylums', 'class']) + '\n')
ofile.write('\t'.join(['FECA-LECA', str(laeca_age), str(leca_age), str(leca.dist), 'Eukaryota', 'Eukaryota', 'Eukaryota']) + '\n')

i = 0
for node in t.traverse():
    if not node.is_leaf():
        death = node.get_distance(node.get_leaf_names()[0])
        length = node.dist
        birth = death + length

        clades = ';'.join(set([x.clade for x in node.get_leaves()]))
        phylum = ';'.join(set([x.phylum for x in node.get_leaves()]))
        classes = ';'.join(set([x.classes for x in node.get_leaves()]))

        if birth >= leca_age and death <= laeca_age and clades == 'Bacteria':
            node.add_feature('name', 'brspace')
            node.add_feature('brspnum', i)

        ofile.write('\t'.join([str(i), str(birth), str(death), str(length), str(clades), phylum, classes]) + '\n')
    else:
        birth = node.dist
        length = node.dist

        clades = ';'.join(set([x.clade for x in node.get_leaves()]))
        phylum = ';'.join(set([x.phylum for x in node.get_leaves()]))
        classes = ';'.join(set([x.classes for x in node.get_leaves()]))

        ofile.write('\t'.join([str(i), str(birth), '0', str(length), str(clades), phylum, classes]) + '\n')

    i += 1
ofile.close()

t.write(features=['name', 'brspnum'], outfile='../outputs/annotated_tree.nwk')