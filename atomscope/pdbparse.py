"""Minimal, dependency-free PDB-format parser.

Reads ATOM/HETATM records from the first model, classifies them into protein
atoms, water, monatomic ions/elements, and organic ligands, and groups hetero
components so they can be analyzed and selected in the UI.
"""

from __future__ import annotations

from dataclasses import dataclass, field

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "SOL"}

STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "MSE", "SEC", "PYL",  # selenomethionine + rare encoded residues
}

# Common biologically relevant metals / monatomic ions seen as HETATM.
METAL_ELEMENTS = {
    "NA", "K", "MG", "CA", "MN", "FE", "CO", "NI", "CU", "ZN", "MO", "W",
    "CD", "HG", "PT", "AU", "AG", "PB", "BA", "SR", "CS", "RB", "LI", "AL",
    "V", "CR", "SN",
}
HALIDE_ELEMENTS = {"CL", "BR", "F", "I"}


@dataclass
class Atom:
    serial: int
    name: str
    res_name: str
    chain: str
    res_seq: int
    icode: str
    x: float
    y: float
    z: float
    element: str
    is_hetero: bool
    occupancy: float = 1.0
    bfactor: float = 0.0

    @property
    def key(self) -> tuple:
        return (self.chain, self.res_seq, self.icode)


@dataclass
class Component:
    """A grouped hetero component (ligand, ion, cofactor, etc.)."""

    res_name: str
    chain: str
    res_seq: int
    icode: str
    atoms: list[Atom] = field(default_factory=list)

    @property
    def label(self) -> str:
        seq = f"{self.res_seq}{self.icode}".strip()
        return f"{self.res_name} {self.chain}{seq}".strip()

    @property
    def kind(self) -> str:
        heavy = [a for a in self.atoms if a.element != "H"]
        if len(heavy) == 1:
            el = heavy[0].element
            if el in METAL_ELEMENTS:
                return "metal"
            if el in HALIDE_ELEMENTS:
                return "ion"
            return "ion"
        return "ligand"


@dataclass
class Structure:
    atoms: list[Atom]
    protein_atoms: list[Atom]
    components: list[Component]
    chains: list[str]

    @property
    def ligand_components(self) -> list[Component]:
        return [c for c in self.components if c.kind == "ligand"]


def _f(s: str, default: float = 0.0) -> float:
    try:
        return float(s)
    except ValueError:
        return default


def parse_pdb(text: str) -> Structure:
    atoms: list[Atom] = []
    seen_alt: dict[tuple, str] = {}
    in_first_model = True

    for line in text.splitlines():
        rec = line[:6].strip()
        if rec == "ENDMDL":
            in_first_model = False
            continue
        if rec == "MODEL":
            # Only keep the first model.
            if atoms:
                in_first_model = False
            continue
        if not in_first_model:
            continue
        if rec not in ("ATOM", "HETATM"):
            continue
        if len(line) < 54:
            continue

        alt_loc = line[16].strip()
        name = line[12:16].strip()
        res_name = line[17:20].strip()
        chain = line[21].strip() or "A"
        res_seq = int(_f(line[22:26], 0))
        icode = line[26].strip()

        # Keep only one alternate location per atom (the first encountered).
        alt_key = (chain, res_seq, icode, name)
        if alt_loc:
            if alt_key in seen_alt and seen_alt[alt_key] != alt_loc:
                continue
            seen_alt[alt_key] = alt_loc

        element = line[76:78].strip().upper() if len(line) >= 78 else ""
        if not element:
            # Fall back: derive from atom name.
            element = "".join(c for c in name if c.isalpha())[:2].upper()

        atoms.append(
            Atom(
                serial=int(_f(line[6:11], 0)),
                name=name,
                res_name=res_name,
                chain=chain,
                res_seq=res_seq,
                icode=icode,
                x=_f(line[30:38]),
                y=_f(line[38:46]),
                z=_f(line[46:54]),
                element=element,
                is_hetero=(rec == "HETATM"),
                occupancy=_f(line[54:60], 1.0) if len(line) >= 60 else 1.0,
                bfactor=_f(line[60:66], 0.0) if len(line) >= 66 else 0.0,
            )
        )

    protein_atoms = [a for a in atoms if not a.is_hetero and a.res_name in STANDARD_AA]

    comp_map: dict[tuple, Component] = {}
    for a in atoms:
        if not a.is_hetero:
            continue
        if a.res_name in WATER_NAMES:
            continue
        key = (a.res_name, a.chain, a.res_seq, a.icode)
        comp = comp_map.get(key)
        if comp is None:
            comp = Component(a.res_name, a.chain, a.res_seq, a.icode)
            comp_map[key] = comp
        comp.atoms.append(a)

    chains = sorted({a.chain for a in protein_atoms})
    components = sorted(
        comp_map.values(), key=lambda c: (c.chain, c.res_seq, c.res_name)
    )
    return Structure(
        atoms=atoms,
        protein_atoms=protein_atoms,
        components=components,
        chains=chains,
    )
