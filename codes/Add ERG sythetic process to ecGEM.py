import cobra
from cobra import Model, Reaction, Metabolite
import os
'''from ETFLdesigner.ETFLdesigner.io.ecModel import load_ecmodel'''
model = cobra.io.read_sbml_model(r'C:\Users\Administrator\Desktop\Modelling\Basic cobra\model\ecYeastGEM_batch.xml')

'''这一步可以开始导入egt模型进入我们的现有模型'''

reaction = Reaction('R_egt11')
reaction.name = 'histidine to Hercynine '
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

SAH_c = Metabolite(
    's_1413[c]',
    formula='C14H20N6O5S',
    name='S-adenosyl-L-homocysteine',
    compartment='c')

Hercynine = Metabolite(
    'Hercynine',
    formula='C9H15N3O2',
    name='Hercynine',
    compartment='c')

prot_egt11 = Metabolite(
    'prot_egt11',
    formula='',
    name='EGT11 enzyme',
    compartment='c')

reaction.add_metabolites({
    model.metabolites.get_by_id('s_1006[c]'): -1.0,  # 组氨酸
    model.metabolites.get_by_id('s_1416[c]'): -3.0,  # SAM
    model.metabolites.get_by_id('s_1413[c]'): 3.0,   # SAH
    Hercynine: 1.0,
    model.metabolites.get_by_id('s_0794[c]'): 3.0,   # H
    prot_egt11: -0.44                    # 这里我们姑且让他和后面的kcat一致，试一试
})
model.add_reactions([reaction])


reaction = Reaction('R_egt12')
reaction.name = 'Hercynine to hercynylcysteine sulfoxide '
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

H2O_c = Metabolite(
    's_0803[c]',
    formula='H2O',
    name='H2O',
    compartment='c')
hercynylcysteine_sulfoxide = Metabolite(
    'hercynylcysteinesulfoxide',
    formula='C12H20N4O5S',
    name='hercynylcysteinesulfoxide',
    compartment='c')
prot_egt12 = Metabolite(
    'prot_egt12',
    formula='',
    name='EGT12 enzyme',
    compartment='c')

reaction.add_metabolites({
    Hercynine: -1.0,
    model.metabolites.get_by_id('s_0981[c]'): -1.0,   # Cys
    model.metabolites.get_by_id('s_1275[c]'): -1.0,    # O2
    model.metabolites.get_by_id('s_0803[c]'): 1.0,    # H2O
    hercynylcysteine_sulfoxide: 1.0,
    prot_egt12: -0.44
})
model.add_reactions([reaction])


reaction = Reaction('R_egt2')
reaction.name = 'hercynylcysteine sulfoxide to Ergothioneine'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = -1000.
reaction.upper_bound = 1000.

erg = Metabolite(
    'ergothioneine',
    formula='C9H15N3O2S',
    name='ergothioneine',
    compartment='c')
prot_egt2 = Metabolite(
    'prot_egt2',
    formula='',
    name='EGT2 enzyme',
    compartment='c')

reaction.add_metabolites({
    model.metabolites.get_by_id('s_0794[c]'): -3.0,   # H2O
    hercynylcysteine_sulfoxide: -1.0,
    model.metabolites.get_by_id('s_1399[c]') : 1.0,   # 丙酮酸
    erg: 1.0,
    prot_egt2: -0.04
})
model.add_reactions([reaction])
model.add_boundary(model.metabolites.get_by_id("ergothioneine"), type="sink")


# Egt2加入蛋白池
reaction = Reaction('draw_EGT2')
reaction.name = 'draw_prot_EGT2'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

prot_egt2 = Metabolite(
    'prot_egt2',
    formula='',
    name='EGT2 enzyme',
    compartment='c')

reaction.add_metabolites({ model.metabolites.get_by_id('prot_pool[c]'): -53.015, prot_egt2: 1})
model.add_reactions([reaction])

# Egt1加入蛋白池
reaction = Reaction('draw_EGT11')
reaction.name = 'draw_prot_EGT11'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ model.metabolites.get_by_id('prot_pool[c]'): -99.025, prot_egt11: 1})
model.add_reactions([reaction])

'''model.remove_reactions(model.reactions.get_by_id('draw_EGT2'))'''


# Egt12加入蛋白池
reaction = Reaction('draw_EGT12')
reaction.name = 'draw_prot_EGT12'
reaction.subsystem = 'Cell cytoplasm Biosynthesis'
reaction.lower_bound = 0.
reaction.upper_bound = 1000.

reaction.add_metabolites({ model.metabolites.get_by_id('prot_pool[c]'): -99.025, prot_egt12: 1})
model.add_reactions([reaction])



output_path = "C:/Users/Administrator/Desktop/ecGEM with ERG.xml"
cobra.io.write_sbml_model(model, output_path)


# constraints添加，以这个为例，是一个减另一个的值为0，即约束二者通量相同
same_flux = model.problem.Constraint(
    model.reactions.____.flux_expression - model.reactions.____.flux_expression, lb=0, ub=0)
model.add_cons_vars(same_flux)

# 对于ecGEM，其constraints直接放在model下，因此直接调整理论上也可以
add_new_into_pool = model.problem.Constraint(model.reactions.prot_pool_exchange.flux_expression -
                                             model.reactions.prot_pool_exchange.flux_expression -
                                             model.reactions.draw_EGT2.flux_expression, lb=0, ub=0)
model.add_cons_vars(add_new_into_pool)




model.remove_reactions(model.reactions.get_by_id('draw_EGT2'))
