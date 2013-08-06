import matplotlib.pyplot as plt
import networkx as nx

proteins = [
	"ACC_pS79",
        "PRAS40_pT246",
	"AKT_pS473",
	"AKT_pT308",
	"BAD_pS112",
	"p27_pT157",
	"p27_pT198",
	"CHK1_pS345",
	"CHK2_pT68",
	"EGFR_pY1068",
	"EGFR_pY1173",
	"EGFR_pY992",
	"4EBP1_pS65",
	"4EBP1_pT37_pT46",
	"HER2_pY1248",
	"HER3_pY1298",
	"ER-alpha_pS118",
	"GSK3-alpha-beta_pS9",
	"GSK3-alpha-beta_pS21_S9",
	"c-JUN_pS73",
	"MEK1_pS217_S221",
	"p38_pT180_Y182",
	"MAPK_pT202_Y204",
	"JNK_pT183_pT185",
	"c-Met_pY1235",
	"mTOR_pS2448",
	"NDRG1_pT346",
	"PDK1_pS241",
	"PEA15_pS116",
	"AMPK_pT172",
	"PKC-alpha_pS657",
	"PKC-pan-betaII_pS660",
	"PKC-delta_pS664",
	"c-Raf_pS338",
	"Rb_pS807_S811",
	"NF-kB-p65_pS536",
	"Rictor_pT1135",
	"p90RSK_pT359_S363",
	"p70S6K_pT389",
	"S6_pS235_S236",
	"S6_pS240_S244",
	"Src_pY416",
	"Src_pY527",
	"STAT3_pY705",
	"YAP_pS127",
	"YB-1_PS102"
        ]

DG = nx.DiGraph()
DG.add_nodes_from(proteins)

#DG.add_edge(, )

## RTK
DG.add_edge("HER2_pY1248", "EGFR_pY1068")
DG.add_edge("HER2_pY1248", "EGFR_pY1173")
DG.add_edge("HER2_pY1248", "EGFR_pY992")
DG.add_edge("HER3_pY1298", "EGFR_pY1068")
DG.add_edge("HER3_pY1298", "EGFR_pY1173")
DG.add_edge("HER3_pY1298", "EGFR_pY992")
DG.add_edge("EGFR_pY1068", "HER2_pY1248")
DG.add_edge("EGFR_pY1068", "HER2_pY1248")
DG.add_edge("EGFR_pY1173", "HER2_pY1248")
DG.add_edge("EGFR_pY1173", "HER3_pY1298")
DG.add_edge("EGFR_pY992", "HER3_pY1298")
DG.add_edge("EGFR_pY992", "HER3_pY1298")
DG.add_edge("EGFR_pY1068", "Src_pY416")
DG.add_edge("EGFR_pY1173", "Src_pY416")
DG.add_edge("EGFR_pY992", "Src_pY416")
DG.add_edge("EGFR_pY1068", "Src_pY527")
DG.add_edge("EGFR_pY1173", "Src_pY527")
DG.add_edge("EGFR_pY992", "Src_pY527")


## MAPK signaling
DG.add_edge("EGFR_pY1068", "c-Raf_pS338")
DG.add_edge("EGFR_pY1173", "c-Raf_pS338")
DG.add_edge("EGFR_pY992", "c-Raf_pS338")
DG.add_edge("c-Met_pY1235", "c-Raf_pS338")
DG.add_edge("EGFR_pY1068", "JNK_pT183_pT185")
DG.add_edge("EGFR_pY1173", "JNK_pT183_pT185")
DG.add_edge("EGFR_pY992", "JNK_pT183_pT185")
DG.add_edge("Src_pY416", "c-Raf_pS338")
DG.add_edge("Src_pY527", "c-Raf_pS338")
DG.add_edge("c-Raf_pS338", "MEK1_pS217_S221")
DG.add_edge("MEK1_pS217_S221", "p38_pT180_Y182")
DG.add_edge("MEK1_pS217_S221", "MAPK_pT202_Y204")
DG.add_edge("JNK_pT183_pT185", "c-JUN_pS73")
DG.add_edge("JNK_pT183_pT185", "PEA15_pS116")
DG.add_edge("MAPK_pT202_Y204", "p90RSK_pT359_S363")
DG.add_edge("MAPK_pT202_Y204", "p70S6K_pT389")
DG.add_edge("p90RSK_pT359_S363", "S6_pS235_S236")
DG.add_edge("p70S6K_pT389", "S6_pS235_S236")
DG.add_edge("p90RSK_pT359_S363", "S6_pS240_S244")
DG.add_edge("p70S6K_pT389", "S6_pS240_S244")


## protein kinase B
DG.add_edge("PDK1_pS241", "PRAS40_pT246")
DG.add_edge("PDK1_pS241", "AKT_pS473")
DG.add_edge("PDK1_pS241", "AKT_pT308")
DG.add_edge("mTOR_pS2448", "PRAS40_pT246")
DG.add_edge("mTOR_pS2448", "AKT_pS473")
DG.add_edge("mTOR_pS2448", "AKT_pT308")
DG.add_edge("Rictor_pT1135", "PRAS40_pT246")
DG.add_edge("Rictor_pT1135", "AKT_pS473")
DG.add_edge("Rictor_pT1135", "AKT_pT308")

DG.add_edge("PRAS40_pT246", "BAD_pS112")

## protein kinase C

## STAT signaling
DG.add_edge("EGFR_pY1068", "STAT3_pY705")
DG.add_edge("EGFR_pY1173", "STAT3_pY705")
DG.add_edge("EGFR_pY992", "STAT3_pY705")
DG.add_edge("HER2_pY1248", "STAT3_pY705")
DG.add_edge("HER3_pY1298", "STAT3_pY705")
DG.add_edge("c-Met_pY1235", "STAT3_pY705")
DG.add_edge("Src_pY416", "STAT3_pY705")
DG.add_edge("Src_pY527", "STAT3_pY705")
DG.add_edge("STAT3_pY705", "NF-kB-p65_pS536")

## hippo signaling

## metabolism/stress 
DG.add_edge("AMPK_pT172", "ACC_pS79")
DG.add_edge("NDRG1_pT346", "NF-kB-p65_pS536")
DG.add_edge("AMPK_pT172", "mTOR_pS2448")

nx.draw(DG)
plt.show()

