"""Protonation and restrained energy minimisation.

Stage 03 places side-chain atoms at ideal internal geometry, so nothing has yet
felt a force field. This module adds hydrogens and relaxes the structure under
a backbone restraint.

What the restraint is for: Route A's central claim is that the peptide backbone
is an experimental observation rather than a model. Restraining backbone heavy
atoms preserves that, while leaving every side chain free to respond to the new
aspartate.

This module does not run dynamics and says nothing about whether the relaxed
conformer is populated in solution. See the stated shortcomings in
docs/design_plan.md.
"""

from __future__ import annotations

from pathlib import Path

import openmm
import openmm.app as app
import openmm.unit as unit
from pdbfixer import PDBFixer

#: amber14 protein parameters with GBn2 implicit solvent. Implicit solvent
#: rather than explicit water because this is a local relaxation, not a
#: simulation, and solvating would multiply the system size for no benefit.
FORCEFIELD_FILES = ("amber14-all.xml", "implicit/gbn2.xml")

#: Backbone heavy atoms. OXT is excluded: it is a terminal carboxylate oxygen
#: whose placement is a modelling choice, not part of the observed trace.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})

#: Harmonic restraint stiffness, kcal/mol/A^2. Firm enough to hold the trace
#: within a few hundredths of an angstrom, loose enough not to fight the
#: force field into numerical trouble.
DEFAULT_RESTRAINT_K = 10.0

#: Non-bonded cutoff for the implicit-solvent system.
NONBONDED_CUTOFF_NM = 1.0


def protonate(path: str | Path, ph: float = 7.4):
    """Add missing heavy atoms and hydrogens to a structure.

    Missing *residues* are deliberately not built. Copy C of 7OW3 begins at
    heavy-chain residue 2, and a design target should not contain invented
    backbone. PDBFixer requires ``missingResidues`` to be set before
    ``addMissingAtoms``, so it is computed and then cleared.

    Returns:
        ``(topology, positions)`` ready for ``ForceField.createSystem``.

    Raises:
        ValueError: if the structure contains a nonstandard residue, which the
            amber14 templates would not cover.
    """
    fixer = PDBFixer(filename=str(path))

    fixer.findMissingResidues()
    fixer.missingResidues = {}

    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        raise ValueError(
            f"{Path(path).name}: nonstandard residues present: {fixer.nonstandardResidues}"
        )

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(ph)
    return fixer.topology, fixer.positions


def backbone_restraint_force(
    topology: app.Topology,
    positions,
    k_kcal_per_mol_per_angstrom2: float = DEFAULT_RESTRAINT_K,
) -> tuple[openmm.CustomExternalForce, list[int]]:
    """Harmonic positional restraint on backbone heavy atoms.

    Returns:
        The force, and the indices of the atoms it restrains.
    """
    force = openmm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    force.addGlobalParameter(
        "k",
        k_kcal_per_mol_per_angstrom2 * unit.kilocalories_per_mole / unit.angstrom**2,
    )
    for name in ("x0", "y0", "z0"):
        force.addPerParticleParameter(name)

    restrained: list[int] = []
    for atom in topology.atoms():
        if atom.name not in BACKBONE_ATOMS:
            continue
        coord = positions[atom.index].value_in_unit(unit.nanometer)
        force.addParticle(atom.index, [coord.x, coord.y, coord.z])
        restrained.append(atom.index)
    return force, restrained


def minimise(
    topology: app.Topology,
    positions,
    restraint_k: float = DEFAULT_RESTRAINT_K,
    platform_name: str = "CPU",
):
    """Minimise with backbone heavy atoms restrained.

    Returns:
        ``(positions, energy_before, energy_after, n_restrained)``. Both
        energies include the restraint term, so they are comparable with each
        other but are not absolute potential energies of the molecule.
    """
    forcefield = app.ForceField(*FORCEFIELD_FILES)
    system = forcefield.createSystem(
        topology,
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=NONBONDED_CUTOFF_NM * unit.nanometer,
        constraints=app.HBonds,
    )
    force, restrained = backbone_restraint_force(topology, positions, restraint_k)
    system.addForce(force)

    # The integrator is required to build a Simulation but is never stepped.
    integrator = openmm.LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = app.Simulation(
        topology, system, integrator, openmm.Platform.getPlatformByName(platform_name)
    )
    simulation.context.setPositions(positions)

    before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    simulation.minimizeEnergy()
    state = simulation.context.getState(getPositions=True, getEnergy=True)
    return (
        state.getPositions(),
        before.value_in_unit(unit.kilocalories_per_mole),
        state.getPotentialEnergy().value_in_unit(unit.kilocalories_per_mole),
        len(restrained),
    )


def write_pdb(topology: app.Topology, positions, path: str | Path) -> None:
    with open(path, "w") as handle:
        app.PDBFile.writeFile(topology, positions, handle, keepIds=True)
