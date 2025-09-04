from __future__ import annotations

CAPABILITIES = {
    "tools": [
        {"name":"featurize_basic",
         "args_schema":{"features":[
             "FG_COUNTS","HALIDE_CLASS","NUCLEOPHILE_CLASS","BORON_CLASS",
             "OXIDATION_DELTA","REDUCTION_DELTA","C_C_FORMATION","C_N_FORMATION",
             "RING_STATS","TPSA","logP"
         ]}},
        {"name":"retrieve_neighbors",
         "args_schema":{"k":"int","filters":{"class_tags":"[string]","condition_core":"[string]"}}},
        {"name":"get_dicts",
         "args_schema":{"reaction_class":"string|null"}},
        {"name":"predict_yield",
         "args_schema":{"reaction":"SMILES","conditions":"object"}},
        {"name":"validate_proposals",
         "args_schema":{"proposals":"object"}}
    ],
    "policies":{"max_actions_per_cycle":3,"max_cycles":4}
}
