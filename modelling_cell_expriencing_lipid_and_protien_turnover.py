
##################################FUNCTIONS#####################################
def modelMatureCell(model,proteinDegrdationRate=1,lipidDegrdationRate=1):

    #constrain light-reactions
    model.reactions.get_by_id("Photon_tx").lower_bound = 0
    model.reactions.get_by_id("Photon_tx").upper_bound = 0

    #constraint sucrose and ammonium uptake
    model.reactions.get_by_id("Sucrose_tx").lower_bound = 0
    model.reactions.get_by_id("Sucrose_tx").upper_bound = 0
    model.reactions.get_by_id("NH4_tx").lower_bound = 0
    model.reactions.get_by_id("NH4_tx").upper_bound = 0

    #set lipid turnover rate
    lipidMW = model.metabolites.get_by_id("PALMITATE_x").formula_weight
    model.reactions.get_by_id("Beta_Oxidation_x").lower_bound = lipidDegrdationRate/lipidMW
    model.reactions.get_by_id("Beta_Oxidation_x").upper_bound = lipidDegrdationRate/lipidMW

    #set protein turnover rate
    proteinMW = model.metabolites.get_by_id("Protein_b").formula_weight
    model.reactions.get_by_id("Protein_degradation").lower_bound = proteinDegrdationRate/proteinMW
    model.reactions.get_by_id("Protein_degradation").upper_bound = proteinDegrdationRate/proteinMW

    #set minimization of glucose uptake as the objective of the  cell
    model.reactions.get_by_id("GLC_tx").objective_coefficient = -1

    #run parsimonious FBA to predict flux distribution with minimal metabolic flux
    sol = flux_analysis.parsimonious.pfba(model)
    model.solution = sol
    return model

#################################### MAIN ##################################
from cobra import io
from cobra import flux_analysis
from cobra.core import Reaction, Metabolite

#import model from SBML file
model = model = io.sbml.create_cobra_model_from_sbml_file("PlantCoreMetabolism_v1_2_3.xml")

# add reaction to represent protein synthesis
rxn = model.reactions.get_by_id("Biomass_tx")
rxn.id = "Protein_biomass"
for met in ["Ca_b","K_b","Mg_b"]:
    met = model.metabolites.get_by_id(met)
    coeff = rxn.metabolites.get(met)
    rxn.add_metabolites({met:-1*coeff})

#create metabolite Protein[b]
met = Metabolite("Protein_b",name="Protein[b]")
formula_dict = rxn.check_mass_balance()
met.formula = "".join([atom+str(formula_dict[atom]*-1) for atom in formula_dict.keys() if atom != "charge"])
met.charge = formula_dict["charge"]*-1

rxn.add_metabolites({met:1})
#check reaction
#print rxn.reaction

# add reaction to represent breakdown of protein to individual amino acids
rxn = Reaction("Protein_degradation")
rxn.add_metabolites({model.metabolites.Protein_b:-1,
                     model.metabolites.L_ALPHA_ALANINE_c:0.2326035201,
                     model.metabolites.ARG_c:0.1771733159,
                     model.metabolites.ASN_c:0.1447158566,
                     model.metabolites.L_ASPARTATE_c:0.3855539754,
                     model.metabolites.GLN_c:0.1325494187,
                     model.metabolites.GLT_c:0.2261465373,
                     model.metabolites.GLY_c:0.2309411553,
                     model.metabolites.HIS_c:0.1673732646,
                     model.metabolites.ILE_c:0.1454255314,
                     model.metabolites.LEU_c:0.0053400685,
                     model.metabolites.LYS_c:0.3415370304,
                     model.metabolites.MET_c:0.1420323181,
                     model.metabolites.PHE_c:0.1980982902,
                     model.metabolites.SER_c:0.3368402626,
                     model.metabolites.THR_c:0.2365218861,
                     model.metabolites.TYR_c:0.2071316851,
                     model.metabolites.VAL_c:0.158571619})
rxn.lower_bound = 0
rxn.upper_bound = 1000
model.add_reaction(rxn)
#check reactions
#print rxn.reaction

#If we assume lipid turnover rate is 0.2 hr−1 and a cell's lipid content is
# 164.80 ng.cell−1, lipid degradation rate can be assumed to be
# 32.96 ng.cell−1.hr−1. Similarly if we assume protein turnover rate is
# 0.5 hr−1 and a cell's protien content is 0.057 ng.cell−1, protein degradation
# rate can be assumed to be 0.029 ng.cell−1.hr
model2 = modelMatureCell(model,proteinDegrdationRate=32.96,lipidDegrdationRate=0.029)

#print summary of flux distribution
model2.summary()

print "From the summary, we can see that cell needs to take up atleast "+ \
        str(model2.reaction.GLC_tx.x)+" nmolglucose.hr−1 \
        and has a respiration rate of "+str(model2.reaction.CO2_tx.x)+ \
        "nmolCO2.hr−1 under optimal flux distribution"

print("Other fluxes of interest:")
rxn = model2.reactions.Mitochondrial_ATP_Synthase_m
print "mitochondrial ATP synthesis flux = "+str(round(rxn.x,2)*3)+" nmolATP/hr"
rxn1 = model2.reactions.GLU6PDEHYDROG_RXN_p
rxn2 = model2.reactions.GLU6PDEHYDROG_RXN_c
print "pentose phosphate pathway flux = "+str(round(rxn1.x+rxn2.x,5))+" nmol G6P/hr"
rxn3 = model2.reactions.PYRUVDEH_RXN_m
print "mitochondrial pyruvate dehydrogenase = "+str(round(rxn3.x,5))+" nmol pyruvate/hr"
rxn4 = model2.reactions.PYRUVDEH_RXN_p
print "plastidial pyruvate dehydrogenase = "+str(round(rxn4.x,5))+" nmol pyruvate/hr"
