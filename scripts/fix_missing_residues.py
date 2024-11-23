#!/usr/bin/env python3

import json
import os
import subprocess
import warnings
import shutil
import string

# Import biopython subpackages
import Bio
import Bio.Align
import Bio.PDB
import Bio.SeqIO


def find_missing(array):
    if array:
        array_set = set(array)
        missing = []
        for i in range(min(array), max(array) + 1):
            if i not in array_set:
                missing.append(i)
        return missing


def parse_structure(pdb_file)
    pdb_id, extension = os.path.splitext(os.path.basename(pdb_file))
    if extension == '.pdb':
        parser = Bio.PDB.PDBParser(QUIET=True)
    if extension == '.cif':
        parser = Bio.PDB.MMCIFParser(QUIET=True)
    reutrn parser.get_structure('structure', pdb_file)


def fix_num_models(pdb_file)
    pdb_data = parse_structure(pdb_file)
    # get the numer of models in the PDB file. 
    num_models = len(pdb_data.child_list)

    # if num_models > 1:
        
    chain_index = 0
    for model in pdb_data:
        for chain in model:
            print(string.ascii_uppercase[chain_index])
                


def check_missing_residues(pdb_file):

    pdb_data = parse_structure(pdb_file)

    has_missing_residues = False
    missing_residues_list = []
    for model in pdb_data:
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

                    # else:
                    #     print(pdb_id, chain.get_id(),
                    #           residue.get_id(), residue.get_resname())

            missing_residues = find_missing(residue_id_list)
            missing_residues_list.append([model.get_id(), chain.get_id(),
                                          missing_residues])

            if missing_residues and len(missing_residues) > 0:
                has_missing_residues = True

    return pdb_id, has_missing_residues, missing_residues_list


def get_sequence(pdb_file, format):
    """Get the amino acid sequence from the the SEQRES records or order
    of atoms in the PDB file """

    if format not in {"pdb-seqres", "pdb-atom"}:
        raise ValueError('format must be pdb-seqres or pdb-atom')
    

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', Bio.BiopythonWarning)
        
        sequence = {}
        for chain in Bio.SeqIO.parse(pdb_file, format):
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

    seqres_sequence = get_sequence(pdb_file, "pdb-seqres")
    atomic_sequence = get_sequence(pdb_file, "pdb-atom")
    
    print(seqres_sequence.keys() == atomic_sequence.keys(), seqres_sequence.keys(), atomic_sequence.keys())

    alignment_list = []
    for chain in sorted(seqres_sequence.keys()):
        target, template = aligner.align(seqres_sequence[chain], atomic_sequence[chain])[0]
        alignment_list.append({
            'target': {'name': 'mytrg', 'seqres': target},
            'template': {'name': f'{name}.{chain}', 'offset': 0, 'seqres': template}
            })

    json_data = {"alignmentlist": alignment_list}
    json.dump(json_data, open(alignment_file, 'w'), indent=4)


def model_structure(pdb_file, output_dir='../output/modeled_structures'):
    """ Model the structure with missing residues using ProMod3
    https://openstructure.org/promod3/3.4/ """

    name = os.path.splitext(os.path.basename(pdb_file))[0]

    alignment_file = os.path.join(output_dir, f'{name}_alignment.json')
    align_pdb_sequences(pdb_file=pdb_file, alignment_file=alignment_file)
    
    shutil.copy(pdb_file, output_dir)
    output_dir = os.path.abspath(output_dir)
    subprocess.run(['docker', 'run', '--rm', '-v', f'{output_dir}:/home',
                    'registry.scicore.unibas.ch/schwede/promod3',
                    'build-model', '-j', f'{name}_alignment.json',
                    '-p', f'{name}.pdb', '-o', f'{name}_model.pdb'
                    ])


if __name__ == '__main__':
    
    import glob
    import pandas as pd
    from tqdm.contrib.concurrent import process_map

    structures = sorted(glob.glob('../data/pdb_downloaded/*.pdb') +
                        glob.glob('../data/pdb_downloaded/*.cif'))

    missing_residues = process_map(check_missing_residues,
                                   structures, chunksize=1)

    df_missing_residues = pd.DataFrame(
        missing_residues, columns=['PDB ID',
                                   'Has Missing Residues',
                                   'Missing Residues'])
    
    # Get the PDB IDs with missing residues
    pdb_ids = df_missing_residues.loc[
        df_missing_residues['Has Missing Residues'], 'PDB ID']
        
    output_directory = '../output/modeled_structures'
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    pdbf_files = [f'../data/pdb_downloaded/{pdb_id}.pdb' for pdb_id in pdb_ids]
    process_map(model_structure, pdbf_files, max_workers=16, chunksize=1)
