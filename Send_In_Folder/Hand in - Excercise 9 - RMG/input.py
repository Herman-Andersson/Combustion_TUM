# To run the script
# conda activate rmg_env
# cd /home/brukare/Combustion/RMG_methane
# python -m rmgpy.rmg.main input.py

# Data sources
database(
    thermoLibraries = ['GRI-Mech3.0', 'primaryThermoLibrary'],
    reactionLibraries = ['GRI-Mech3.0',"Klippenstein_Glarborg2016"],
    seedMechanisms = ['GRI-Mech3.0'],
    kineticsDepositories = ['training'],
    kineticsFamilies = 'default',
    kineticsEstimator = 'rate rules',
)

# List of species
species(
    label='CH4',
    reactive=True,
    structure=SMILES("C"),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)

# Reactor conditions to ensure comprehensive mechanism
simpleReactor(
    temperature=[(700,'K'), (2000,'K')],  # Only 2 temperatures
    pressure=(1,'bar'),
    nSims=10,
    initialMoleFractions={
        "CH4": [0.0499, 0.1361],  # φ = 0.5 and 1.5 
        "O2": [0.1815, 0.1996],   # Corresponding O2
        "N2": 0.71645    # Corresponding N2
    },
    terminationConversion={'CH4': 0.99},
    terminationTime=(1e6,'s'),
    balanceSpecies = "N2",
)


# Simulator tolerances
simulator(
    atol=1e-16,
    rtol=1e-8,
)

# Model tolerances
model(
    toleranceMoveToCore=0.1,
    toleranceInterruptSimulation=0.1,
    maximumEdgeSpecies=5000,
)

# Pressure dependence
pressureDependence(
    method='modified strong collision',
    maximumGrainSize=(0.5,'kcal/mol'),
    minimumNumberOfGrains=250,
    temperatures=(300,3000,'K',8),
    pressures=(0.001,100,'bar',5),
    interpolation=('PDepArrhenius', 6, 4),
)

# Constraints on generated species
generatedSpeciesConstraints(
    maximumCarbonAtoms = 3,
    allowed=['input species','seed mechanisms','reaction libraries'],
)

# Miscellaneous options
options(
    units='si',
    generateOutputHTML=True,
    generatePlots=False,
    saveSimulationProfiles=True,
    verboseComments=False,
    saveEdgeSpecies=True,
)