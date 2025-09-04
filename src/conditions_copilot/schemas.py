from __future__ import annotations
from typing import List, Dict, Any, Optional
from pydantic import BaseModel, Field, ConfigDict

class Neighbor(BaseModel):
    ConditionCore: str
    tails: str
    yield_: float = Field(..., alias="yield")
    substrate_tags: List[str]
    notes: Optional[str] = None

class Coverage(BaseModel):
    n_total: int
    n_close: int
    chemotype_overlap: str
    median_yield_close: Optional[float] = None

class Retrieval(BaseModel):
    neighbors: List[Neighbor]
    coverage: Coverage
    negatives: List[Dict[str, Any]] = []
    note: Optional[str] = None

class Dictionaries(BaseModel):
    catalysts: List[str]
    ligands: List[str]
    bases: List[str]
    solvents: List[str]
    additives: List[str]

class Constraints(BaseModel):
    metals_allowed: List[str] = ["Pd","Ni","Cu"]
    metal_blacklist: List[str] = []
    solvent_blacklist: List[str] = []
    env: str = "inert_required"
    temp_limits_C: Dict[str,int] = {"min":20,"max":120}
    time_limits_h: Dict[str,int] = {"max":20}
    safety_flags: List[str] = ["no_DMF_at>120C"]

class Features(BaseModel):
    model_config = ConfigDict(extra="ignore")
    HALIDE_CLASS: Optional[str] = None
    NUCLEOPHILE_CLASS: Optional[str] = None
    BORON_CLASS: Optional[str] = None
    OXIDATION_DELTA: Optional[bool] = None
    REDUCTION_DELTA: Optional[bool] = None
    C_C_FORMATION: Optional[bool] = None
    C_N_FORMATION: Optional[bool] = None
    ORTHO_COUNT: Optional[int] = None
    EWG_FLAGS: Optional[List[str]] = None
    RING_STATS: Optional[Dict[str,int]] = None
    TPSA: Optional[float] = None
    logP: Optional[float] = None

class Reaction(BaseModel):
    smiles: str
    type_hint: Optional[str] = None
    class_confidence: Optional[float] = None

class DiscoveryPayload(BaseModel):
    reaction: Reaction
    context: Dict[str, Any]
    capabilities: Dict[str, Any]
    constraints: Constraints
    state: Dict[str, Any]
    ask: str
    dictionaries: Optional[Dictionaries] = None

class ProposalsRequest(BaseModel):
    reaction: Reaction
    features: Features
    retrieval: Retrieval
    dictionaries: Dictionaries
    constraints: Constraints
    ask: str

class Reagent(BaseModel):
    role: str
    name: str
    amount_mol_: Optional[float] = Field(None, alias="amount_mol%")
    amount_equiv: Optional[float] = None
    volume_mL_per_mmol: Optional[float] = None

class FullConditions(BaseModel):
    reagent_order: List[Reagent]
    temperature_C: Dict[str, Any]
    time_h: Dict[str, Any]
    concentration_M: Optional[Dict[str, Any]] = None
    atmosphere: Optional[str] = None
    notes: Optional[List[str]] = None

class Proposal(BaseModel):
    condition_core: str
    full_conditions: FullConditions
    applicability_tags: List[str]
    support: Dict[str, Any]
    risks: List[str]
    confidence: float

class ScreeningPanel(BaseModel):
    when_to_use: str
    design: Dict[str, Any]

class ProposalsResponse(BaseModel):
    proposals: List[Proposal]
    screening_panel: Optional[ScreeningPanel] = None
    validation: Optional[Dict[str, Any]] = None
    next_actions: Optional[List[Dict[str, Any]]] = None
