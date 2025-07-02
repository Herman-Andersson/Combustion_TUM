import os
import cantera as ct

# Define a gas mixture at a high temperature that will undergo a reaction:
gas = ct.Solution('gri30.yaml')
# Set equivalence ratio to 1 for stoichiometric methane-air combustion
gas.set_equivalence_ratio(1.0, 'CH4', 'O2:2.0, N2:7.52')
gas.TP = 1500, 1e5  # 1500 K, 1 bar (100000 Pa)

# Define the element to follow in the reaction path diagram:
element = 'C'

# Define a reactor, let it react until the temperature reaches 1800 K:
r = ct.IdealGasReactor(gas)
net = ct.ReactorNet([r])
T = r.T
while T < 1800:
    net.step()
    T = r.T

# Initiate the reaction path diagram:
diagram = ct.ReactionPathDiagram(gas, element)

# Options for cantera:
diagram.show_details = False
diagram.font='CMU Serif Roman'
diagram.threshold=0.02 # Thresshold for minimum flux relative value
diagram.dot_options='node[fontsize=20,shape="box"]'
diagram.title = 'Reaction path diagram following {0}'.format(element)


# Define the filenames:
dot_file = 'ReactionPathDiagram.dot'
img_file = 'ReactionPathDiagram.png'

# Write the dot-file first, then create the image from the dot-file with customizable
# parameters:
diagram.write_dot(dot_file)

# The command -Tpng defines the filetype and needs to fit your filename defined above,
# or else you will get errors opening the file later.
# The command -Gdpi sets the resolution of the generated image in dpi.
# os.system('dot {0} -Tpng -o{1} -Gdpi=300'.format(dot_file, img_file))
os.system('dot {0} -Earrowhead="onormal" -Esamehead="true" -Nstyle="filled" -Nfillcolor="lightgreen" -Nshape="box3d" -Tpng -Gdpi=300 -Gbgcolor="#FFFFFF" -Gnodesep=0.1 -o{1}'.format(dot_file, img_file))


print(f"Reaction path diagram saved as: {img_file}")
print("Done!")