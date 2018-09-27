# requires trimesh (install using "conda install -c conda-forge trimesh")
import numpy as np
import trimesh
from mercury_interface import MercuryInterface
helper = MercuryInterface()
entry = helper.current_entry
mol = entry.molecule
atom_array = []
atom_type_set = set()
stl_dict = {}

dummy = trimesh.primitives.Sphere(radius=0.0, center=(0.,0.,0.))
SCALE_FACTOR = 10.0  # scale factor in mm/Angstrom
html_file_name = helper.output_html_file
html_file = open(html_file_name, "w")

for atom in mol.atoms:
    atom_type_set.add(atom.atomic_symbol)

stl_dict = dict.fromkeys(atom_type_set, dummy)

for atom in mol.atoms:
    stl_dict[atom.atomic_symbol] += trimesh.primitives.Sphere(radius=atom.vdw_radius*0.25*SCALE_FACTOR, center=np.multiply(atom.coordinates, SCALE_FACTOR))
    for bond in atom.bonds:
        atom1, atom2 = bond.atoms
        if atom1.label == atom.label:               # always set p1 to atom under consideration
            p1 = np.multiply(np.array(atom1.coordinates), SCALE_FACTOR)
            p2 = np.multiply(np.array(atom2.coordinates), SCALE_FACTOR)
        elif atom2.label == atom.label:
            p1 = np.multiply(np.array(atom2.coordinates), SCALE_FACTOR)
            p2 = np.multiply(np.array(atom1.coordinates), SCALE_FACTOR)
        mp = (p1 + p2)/2.0                          # midpoint of bond
        z = np.array([0., 0., 1.])                  # unit vector along z (initial direction of cylinder)
        d1 = np.subtract(p1, mp)                    # difference between atom and midpoint
        length = np.linalg.norm(d1)                 # distance between atom and midpoint
        t = np.cross(z, d1)                         # t-vector normal to bond/z plane
        angle = np.arccos(np.dot(z, d1) / length)   # angle to rotate bond to correct orientation
        rot = trimesh.transformations.rotation_matrix(angle, t)    # rotation transformation
        trans1 = trimesh.transformations.translation_matrix(-d1/2.0+p1)  # translation transformation
        tf1 = trans1.dot(rot)                       # combined transformation
        stl_dict[atom.atomic_symbol] += trimesh.primitives.Cylinder(radius=0.20*SCALE_FACTOR, height=length, sections=20, transform=tf1)

for key,value in stl_dict.items():
    html_file.write('Writing mesh file for ' + str(key) + '<br />')
    stl_file_name = str(entry.identifier) + '_' + str(key) + '.stl'
    stl_file = open(stl_file_name, 'w')
    stl_mesh = trimesh.io.export.export_stl_ascii(value)
    stl_file.write(stl_mesh)
    stl_file.close()


html_file.write('Mesh generation complete<br />')
html_file.close()
