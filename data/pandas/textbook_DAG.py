import matplotlib.pyplot as plt
import networkx as nx

proteins = [
	"ACACA_pS79.ACACB_pS79",
        "AKT1S1_pT246",
	"AKT1_pS473.AKT2_pS473.AKT3_pS473",
	"AKT1_pT308.AKT2_pT308.AKT3_pT308",
	"BAD_pS112",
	"CDKN1B_pT157",
	"CDKN1B_pT198",
	"CHEK1_pS345",
	"CHEK2_pT68",
	"EGFR_pY1068",
	"EGFR_pY1173",
	"EGFR_pY992",
	"EIF4EBP1_pS65",
	"EIF4EBP1_pT37_pT46",
	"ERBB2_pY1248",
	"ERBB3_pY1298",
	"ESR1_pS118",
	"GSK3A_pS9.GSK3B_pS9",
	"GSK3A_pS9_pS21.GSK3B_pS9_pS21",
	"JUN_pS73",
	"MAP2K1_pS217_S221",
	"MAPK14_pT180_Y182",
	"MAPK1_pT202_Y204.MAPK3_pT202_Y204",
	"MAPK8_pT183_pT185",
	"MET_pY1235",
	"MTOR_pS2448",
	"NDRG1_pT346",
	"PDK1_pS241",
	"PEA15_pS116",
	"PRKAA1_pT172.PRKAA2_pT172",
	"PRKCA_pS657",
	"PRKCA_pS660.PRKCB_pS660.PRKCD_pS660.PRKCE_pS660.PRKCH_pS660.PRKCQ_pS660",
	"PRKCD",
	"RAF1_pS338",
	"RB1_pS807_S811",
	"RELA_p65_pS536",
	"RICTOR_pT1135",
	"RPS6KA1_pT359_S363",
	"RPS6KB1_pT389",
	"RPS6_pS235_S236",
	"RPS6_pS240_S244",
	"SRC_pY416",
	"SRC_pY527",
	"STAT3_pY705",
	"WWTR1_pS89",
	"YAP1_pS127",
	"YBX1_PS102"
        ]

DG = nx.DiGraph()
DG.add_nodes_from(proteins)

#DG.add_edge(, )

## RTK
DG.add_edge("EGFR_pY1068", "SRC_pY416")
DG.add_edge("EGFR_pY1173", "SRC_pY416")
DG.add_edge("EGFR_pY992", "SRC_pY416")
DG.add_edge("ERBB2_pY1248", "SRC_pY416")
DG.add_edge("ERBB3_pY1298", "SRC_pY416")
DG.add_edge("EGFR_pY1068", "SRC_pY527")
DG.add_edge("EGFR_pY1173", "SRC_pY527")
DG.add_edge("EGFR_pY992", "SRC_pY527")
DG.add_edge("ERBB2_pY1248", "SRC_pY527")
DG.add_edge("ERBB3_pY1298", "SRC_pY527")


## MAPK signaling
DG.add_edge("EGFR_pY1068", "RAF1_pS338")
DG.add_edge("EGFR_pY1173", "RAF1_pS338")
DG.add_edge("EGFR_pY992", "RAF1_pS338")
DG.add_edge("ERBB2_pY1248", "RAF1_pS338")
DG.add_edge("ERBB3_pY1298", "RAF1_pS338")
DG.add_edge("MET_pY1235", "RAF1_pS338")
DG.add_edge("SRC_pY416", "RAF1_pS338")
DG.add_edge("SRC_pY527", "RAF1_pS338")
DG.add_edge("RAF1_pS338", "MAP2K1_pS217_S221")
DG.add_edge("MAP2K1_pS217_S221", "MAPK14_pT180_Y182")
DG.add_edge("MAP2K1_pS217_S221", "MAPK1_pT202_Y204.MAPK3_pT202_Y204")
DG.add_edge("MAP2K1_pS217_S221", "MAPK8_pT183_pT185")
DG.add_edge("MAPK8_pT183_pT185", "JUN_pS73")
DG.add_edge("MAPK1_pT202_Y204.MAPK3_pT202_Y204", "RPS6KA1_pT359_S363")
DG.add_edge("MAPK1_pT202_Y204.MAPK3_pT202_Y204", "RPS6KB1_pT389")
DG.add_edge("RPS6KA1_pT359_S363", "RPS6_pS235_S236")
DG.add_edge("RPS6KB1_pT389", "RPS6_pS235_S236")
DG.add_edge("RPS6KA1_pT359_S363", "RPS6_pS240_S244")
DG.add_edge("RPS6KB1_pT389", "RPS6_pS240_S244")


## protein kinase B
DG.add_edge("PDK1_pS241", "AKT1S1_pT246")
DG.add_edge("PDK1_pS241", "AKT1_pS473.AKT2_pS473.AKT3_pS473")
DG.add_edge("PDK1_pS241", "AKT1_pT308.AKT2_pT308.AKT3_pT308")
DG.add_edge("MTOR_pS2448", "AKT1S1_pT246")
DG.add_edge("MTOR_pS2448", "AKT1_pS473.AKT2_pS473.AKT3_pS473")
DG.add_edge("MTOR_pS2448", "AKT1_pT308.AKT2_pT308.AKT3_pT308")
DG.add_edge("RICTOR_pT1135", "AKT1S1_pT246")
DG.add_edge("RICTOR_pT1135", "AKT1_pS473.AKT2_pS473.AKT3_pS473")
DG.add_edge("RICTOR_pT1135", "AKT1_pT308.AKT2_pT308.AKT3_pT308")

DG.add_edge("AKT1S1_pT246", "BAD_pS112")

## protein kinase C

## STAT signaling
DG.add_edge("EGFR_pY1068", "STAT3_pY705")
DG.add_edge("EGFR_pY1173", "STAT3_pY705")
DG.add_edge("EGFR_pY992", "STAT3_pY705")
DG.add_edge("ERBB2_pY1248", "STAT3_pY705")
DG.add_edge("ERBB3_pY1298", "STAT3_pY705")
DG.add_edge("MET_pY1235", "STAT3_pY705")
DG.add_edge("SRC_pY416", "STAT3_pY705")
DG.add_edge("SRC_pY527", "STAT3_pY705")


## hippo signaling
DG.add_edge("YAP1_pS127", "WWTR1_pS89")

## metabolism 
DG.add_edge("PRKAA1_pT172.PRKAA2_pT172", "ACACA_pS79.ACACB_pS79")


nx.draw(DG)
plt.show()

