#!/usr/bin/env python3

import collections
import glob
import json
import os
import shutil
import string
import subprocess
import warnings

# Import biopython subpackages
import Bio
import Bio.Align
import Bio.PDB
import Bio.SeqIO
import pandas as pd
from tqdm.contrib.concurrent import process_map


def make_directory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def find_missing(array):
    if array:
        array_set = set(array)
        missing = []
        for i in range(min(array), max(array) + 1):
            if i not in array_set:
                missing.append(i)
        return missing


def parse_structure(pdb_file):
    pdb_id, extension = os.path.splitext(os.path.basename(pdb_file))
    if extension == '.pdb':
        parser = Bio.PDB.PDBParser(QUIET=True)
    if extension == '.cif':
        parser = Bio.PDB.MMCIFParser(QUIET=True)
    return parser.get_structure('pdb_id', pdb_file)


def merge_models(pdb_file, output_file):
    """ Provide structure file path """
    # Read the PDB or mmCIF file
    structure = parse_structure(pdb_file)
    
    seqres_records = collections.defaultdict(list)
    with open(pdb_file, 'r') as in_file:
        for line in in_file:
            if line.startswith('SEQRES'):
                chain = line[11]
                seqres_records[chain].append(line)

    seqres_output = ''
    index = 0
    for model in structure:
        for chain in model:
            old_chain_name = chain.get_id()
            new_chain_name = string.ascii_uppercase[index]
            chain.id = new_chain_name

            # Add a SEQRES record for the chain
            for seqres in seqres_records[old_chain_name]:
                seqres_output += seqres[:11] + new_chain_name + seqres[12:]

            if model.get_id() != 0:
                chain.detach_parent()
                structure[0].add(chain)
            index += 1

    # Write the modified PDB file
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_file)
    
    # Add the SEQRES records to the modified PDB file
    with open(output_file, 'r+') as out_file:
        content = out_file.read()
        out_file.seek(0, 0)
        out_file.write(seqres_output.rstrip('\r\n') + '\n' + content)


def check_missing_residues(pdb_file):
    """ Provide biopython PDB structure Bio.PDB.PDBParser or Bio.PDB.MMCIFParser
    """
    # Read the PDB or mmCIF file
    structure = parse_structure(pdb_file)

    has_missing_residues = False
    list_missing_residues = []
    for model in structure:
        for chain in model:
            residue_id_list = []
            for residue in chain:
                if residue.get_id()[0] == ' ':  # Ignore water and heteroatoms
                    if residue.get_id()[2] == ' ':  # Check insertion code
                        atom_names = {atom.get_name() for atom in residue}

                        # # Check if all backbone atoms are present in the residue
                        # if {'N', 'CA', 'C', 'O'}.issubset(atom_names):

                        # Check at least C-alpha atom for the residue is present
                        if {'CA'}.issubset(atom_names):
                            residue_id_list.append(residue.get_id()[1])
            missing_residues = find_missing(residue_id_list)
            list_missing_residues.append([model.get_id(), chain.get_id(),
                                          missing_residues])

            if missing_residues and len(missing_residues) > 0:
                has_missing_residues = True

    return has_missing_residues, list_missing_residues
    

def get_sequence(pdb_file, format):
    """Get the amino acid sequence from the the SEQRES records or order
    of atoms in the PDB file """

    name, ext = os.path.splitext(os.path.basename(pdb_file))
    
    if format not in {"seqres", "atom"}:
        raise ValueError("format must be 'seqres' or 'atom'")

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', Bio.BiopythonWarning)
        
        sequence = {}
        for chain in Bio.SeqIO.parse(pdb_file, format=f'{ext[1:]}-{format}'):
            sequence[chain.annotations['chain']] = chain.seq.replace('X', '')
        return sequence


def align_pdb_sequences(pdb_file, alignment_file):
    """Align the amino acid sequence from the the SEQRES records to the
    sequence from the order of atoms in the PDB file."""

    name = os.path.splitext(os.path.basename(pdb_file))[0]

    # Calling the biopython alignment method
    aligner = Bio.Align.PairwiseAligner(
        open_gap_score = -0.5,
        extend_gap_score = -0.1,
        target_end_gap_score = 0.0,
        query_end_gap_score = 0.0
        )

    seqres_sequence = get_sequence(pdb_file, "seqres")
    atomic_sequence = get_sequence(pdb_file, "atom")
    
    print(name, seqres_sequence.keys() == atomic_sequence.keys(), seqres_sequence.keys(), atomic_sequence.keys())

    alignment_list = []
    for chain in sorted(atomic_sequence.keys()):
        target, template = aligner.align(seqres_sequence[chain], atomic_sequence[chain])[0]
        alignment_list.append({
            'target': {'name': 'mytrg', 'seqres': target},
            'template': {'name': f'{name}.{chain}', 'offset': 0, 'seqres': template}
            })

    json_data = {"alignmentlist": alignment_list}
    json.dump(json_data, open(alignment_file, 'w'), indent=4)


def model_structure(pdb_file, output_dir='./'):
    """ Model the structure with missing residues using ProMod3
    https://openstructure.org/promod3/3.4/ """

    name = os.path.splitext(os.path.basename(pdb_file))[0]
    
    alignment_file = os.path.join(output_dir, f'{name}_alignment.json')
    align_pdb_sequences(pdb_file=pdb_file, alignment_file=alignment_file)

    shutil.copy(pdb_file, output_dir)

    template_structure = os.path.basename(pdb_file)
    output_path = os.path.abspath(output_dir)
    subprocess.run(['docker', 'run', '--rm', '-v', f'{output_path}:/home',
                    'registry.scicore.unibas.ch/schwede/promod3',
                    'build-model', '-j', f'{name}_alignment.json',
                    '-p', template_structure, '-o', f'{name}_model.pdb'
                    ])


def main(pdb_file):
    pdb_id = os.path.splitext(os.path.basename(pdb_file))[0]
    
    if os.path.exists(f'../output/modeled_structures/{pdb_id}_model.pdb'):
        return
    else:
        try:
            # Check if multiple models are present in the PDB file.
            structure = parse_structure(pdb_file)

            has_missing_res, _ = check_missing_residues(pdb_file)
            
            num_models = len(structure.child_list)
            if num_models > 1:
                input_pdb_file = os.path.join(f'../output/merged_structures/{pdb_id}.pdb')
                merge_models(
                    pdb_file, 
                    output_file=input_pdb_file
                )
            else:
                input_pdb_file = pdb_file

            if has_missing_res:
                model_structure(input_pdb_file, output_dir='../output/modeled_structures')
        except:
            with open('../output/failed_pdbs.txt', 'a') as failed_pdbs:
                print(pdb_id, file=failed_pdbs)

if __name__ == '__main__':
    pdb_structures = sorted(glob.glob('../data/pdb_downloaded/*.pdb')) # +
                            # glob.glob('../data/pdb_downloaded/*.cif'))

    make_directory('../output/modeled_structures')
    make_directory('../output/merged_structures')

    for pdb_structure in pdb_structures:
        main(pdb_structure)
    
    # main('../data/pdb_downloaded/3bcr.cif')
